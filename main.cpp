
/**
 * MALVA - genotyping by Mapping-free ALternate-allele detection of known
 *VAriants Copyright (C) 2019  Giulia Bernardini, Luca Denti, Marco Previtali
 *
 * This file is part of MALVA.
 *
 * MALVA is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MALVA is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MALVA; see the file LICENSE. If not, see
 * <https://www.gnu.org/licenses/>.
 **/

#include <algorithm>
//#include <chrono>
//#include <fstream>
#include <iostream>
#include <utility>

#include <set>
#include <string>
#include <unordered_map>
#include <vector>
#include <string_view>

#include <math.h>
#include <zlib.h>

#include <omp.h>
//#include <openacc.h>

#include <future>
#include <utility>

#include "hts_log.h"
#include "kmc_api/kmc_file.h"
#include "kseq.h"
#include "vcf.h"

#include "argument_parser.hpp"
#include "bloom_filter.hpp"
#include "kmap.hpp"
#include "var_block.hpp"


#include "MurmurHash3.hpp"


void pelapsed(const std::string &s , double start_t) {
    auto now_t = omp_get_wtime();
    std::cerr << "[malva-geno/" << s << "] Time elapsed "
              << now_t -start_t << "s" << std::endl;
}

KSEQ_INIT(gzFile, gzread)

extern const char RCN[128];

//#pragma acc declare create(RCN)

std::unordered_map<std::string, std::string> read_references()
{
    gzFile fasta_in = gzopen(opt::fasta_path.c_str(), "r");
    kseq_t *reference = kseq_init(fasta_in);
    std::unordered_map<std::string, std::string> refs;
    while (kseq_read(reference) >= 0)
    {
        std::string id = reference->name.s;
        if (id.compare(0, 3, "chr") == 0) {
            id = id.substr(3);
        }
        //std::string seq(reference->seq.s);
        //std::transform(std::execution::par_unseq, seq.begin(), seq.end(), seq.begin(), ::toupper);
        refs[id] = reference->seq.s;
    }

    kseq_destroy(reference);
    gzclose(fasta_in);
    return refs;
}

void upper_sequence(std::string& seq)
{
    #pragma omp parallel for simd
    for(long unsigned int i = 0; i < seq.size(); ++i)
        seq[i] &= ~0x20;
}
std::vector<Variant> read_variants()
{
    htsFile *vcf = bcf_open(opt::vcf_path.c_str(), "r");
    bcf_hdr_t *vcf_header = bcf_hdr_read(vcf);
    int is_file_flag = 0;
    if(opt::samples != "-")
        is_file_flag = 1;
    int set_samples_code = bcf_hdr_set_samples(vcf_header, opt::samples.c_str(), is_file_flag);
    if(set_samples_code != 0) {
        throw "ERROR: VCF samples subset invalid code\n";
    }
    bcf1_t *vcf_record = bcf_init();
    std::vector<Variant> vs;
    while (bcf_read(vcf, vcf_header, vcf_record) == 0)
    {
        bcf_unpack(vcf_record, BCF_UN_STR);
        vs.emplace_back(vcf_header, vcf_record, opt::freq_key);
    }
    bcf_hdr_destroy(vcf_header);
    bcf_destroy(vcf_record);
    bcf_close(vcf);
    return vs;
}
/**
 * Method to add kmers to the bloom filter
 **/
void add_kmers_to_bf(BF &bf, KMAP &ref_bf, const VK_GROUP &kmers)
{
    //omp_set_nested(0);
    //#pragma omp parallel if(omp_get_active_level() < 1)
    {
        //#pragma omp single nowait
        {
            for (const auto& v : kmers)
                // For each variant
                for (const auto &p : v.second)
                    // For each allele of the variant/
                    for (const auto &Ks : p.second)
                        // For each list of kmers of the allele
                        for (const std::string& kmer : Ks)
                        {
                            // For each kmer in the kmer list
                            //#pragma omp task //if(omp_in_parallel())
                            {
                                if (p.first == 0)
                                    ref_bf.add_key(kmer);
                                else
                                    bf.add_key(kmer);
                            }
                        }
            }
        }
    }
// potential data race/memory leak valgrind??
void compute_main_filters(std::vector<Variant> &vs, BF &bf, KMAP &ref_bf,
                       const std::unordered_map<std::string, std::string> &refs)
{
    std::string last_seq_name = "";
    VB vb(opt::k, opt::error_rate);
    #pragma omp parallel //proc_bind(close)
    {
        #pragma omp single 
        {
            for (long unsigned int i = 0; i < vs.size(); ++i)
            {
                Variant v = vs.at(i);
                if (last_seq_name.size() == 0)
                    last_seq_name = v.seq_name;
                // for(auto v = vs.begin(); v != vs.end() ; ++v){
                if (!v.has_alts or !v.is_present)
                    continue;
                if (vb.empty())
                {
                    vb.add_variant(v);
                    continue;
                }
                if (!vb.is_near_to_last(v) || last_seq_name != v.seq_name)
                {
                    //std::cerr<<"size notask "<<vb.get_size()<<std::endl;
                    #pragma omp task firstprivate(vb, last_seq_name) untied
                    { 
                        std::string ref = refs.at(last_seq_name);
                        VK_GROUP kmers = vb.extract_kmers(ref);
                        add_kmers_to_bf(bf, ref_bf, kmers);
                    }
                    vb.clear();
                    if (last_seq_name != v.seq_name)
                        last_seq_name = v.seq_name;
                }
                vb.add_variant(v);
            }
            if (!vb.empty())
            {
                VK_GROUP kmers = vb.extract_kmers(refs.at(last_seq_name));
                add_kmers_to_bf(bf, ref_bf, kmers);
            }
        }
    }
}

void compute_context_filter(const BF &bf, BF &context_bf,
                         const std::unordered_map<std::string, std::string> &refs,
                         const std::vector<std::string> &used_seq_names) {
    omp_set_nested(1);
    if(used_seq_names.size() == 0)
        return;
    for (long unsigned int seq = 0; seq < used_seq_names.size(); ++seq)
    {
        long unsigned int size = refs.at(used_seq_names[seq]).size();
        const char* ref = refs.at(used_seq_names[seq]).c_str();
        int pos = (opt::ref_k - opt::k) / 2;
        //#pragma omp parallel for proc_bind(spread)#pragma omp end declare target
        //std::cerr<<"num devices: "<<omp_get_num_devices()<<std::endl;
        uint64_t* cv = new uint64_t[size];
        //#pragma acc parallel loop  independent gang copyin(ref[:size], pos) copyout(cv[:size])
        #pragma omp target teams distribute parallel for map(to: ref[0:size], pos, RCN[0:128]) map(cv[0:size])
        for (long unsigned int p = opt::ref_k; p < size; ++p)
        {
            auto it0 = p - opt::ref_k;
            auto it1 = p - opt::ref_k + pos;
            char t[64];
            char ckmer[64];
            char* canon;
            #pragma omp simd
            for(int i =0 ; i < opt::k ; ++i)
                t[i] = ref[it1 + i];
            //cv[it0]=calc_hash(t, opt::k);
            //reverse_cmpl(ckmer, kmer, k);
            #pragma omp simd
            for(int i = 0; i < opt::k; ++i) {
                ckmer[i] = RCN[(int)t[opt::k - 1 - i]];
            }
            
            int i = 0;
            //#pragma acc loop seq private(i)
            while(i < opt::k && t[i] == ckmer[i])
                ++i;
            
            uint64_t hash[2];
            if(i == opt::k || t[i] < ckmer[i])
                //randomfun(5);
                canon = t;
            else
                canon = ckmer;
            MurmurHash3_x64_128(canon, opt::k, 0, reinterpret_cast<void *>(&hash));
                //  else
            //randomfun(canon, 6);
                //MurmurHash3_x64_128(ckmer, opt::k, 0, reinterpret_cast<void *>(&hash));
            cv[it0] = hash[0];	
            //if(bf.test_key(std::string_view(&reference[it1], opt::k)))
            //  context_bf.add_key(std::string_view(&reference[it0], opt::ref_k));
            
        }
        #pragma omp parallel for
        for(long unsigned int i = 0; i < size; ++i)
        {
            if(bf.ttest(cv[i]))
                context_bf.add_key(std::string_view(&ref[i], opt::ref_k));
        }
        delete[] cv;
    }
}
//int d = ceil(size/(1<<20));
            /*for(int i = 0 ; i < d ; ++i)
        {
            long unsigned int lb = opt::ref_k + ((1<<20) * i);
            long unsigned int ub = (opt::ref_k + ((1<<20) * (i+1)));
            if(ub > size)
                ub =size;
            */


/**
 * Method to compute and store the coverages of the alleles of the
 * variants of a var_block. It uses the coverages stored in the bloom
 * filters/map.
 *	std::cerr<<"omp sections done "<<t2-i*/
void set_coverages(BF &bf, KMAP &ref_bf, VB &vb,
                   const VK_GROUP &kmers /*, const float &cap*/)
{
    for (const auto &var : kmers)
    {
        // For each variant
        Variant v = vb.get_variant(var.first);
        for (const auto &p : var.second)
        {
            float allele_cov = 0;
            for (const auto &Ks : p.second)
            {
                float curr_cov = 0;
                int n = 0; // Number of kmers in the signature
                for (const auto &kmer : Ks)
                {
                    int w = 0;
                    if (p.first == 0)
                        w = ref_bf.get_count(kmer.c_str());
                    else
                        w = bf.get_count(kmer.c_str());
                    if (w > 0) { // maybe useless
                        curr_cov = (curr_cov * n + w) / (n + 1);
                        ++n;
                    }
                }
                if (curr_cov > allele_cov)
                    allele_cov = curr_cov;
            }
            // we can now set the allele coverage
            vb.set_variant_coverage(var.first, p.first, allele_cov);
        }
    }
}

/**
 * Method to clean and print VCF header. It adds GT and GQ FORMAT,
 * removes all samples, and adds donor sample.
 **/
void print_cleaned_header()
{
    htsFile *vcf = bcf_open(opt::vcf_path.c_str(), "r");
    bcf_hdr_t *vcf_header = bcf_hdr_read(vcf);
    bcf1_t *vcf_record = bcf_init();

    // Adding format fields - if already present, they won't be added
    bcf_hdr_append(vcf_header, "##FORMAT=<ID=GT,Number=1,Type=String,"
                   "Description=\"Genotype\">");
    bcf_hdr_append(vcf_header, "##FORMAT=<ID=GQ,Number=1,Type=Integer,"
                   "Description=\"Genotype Quality\">");

    // Adding donor sample and removing all other samples
    const char *new_sample = "DONOR";
    bcf_hdr_add_sample(vcf_header, new_sample);
    bcf_hdr_sync(vcf_header);
    bcf_hdr_set_samples(vcf_header, new_sample, 0);

    // Formatting and printing header
    kstring_t htxt = {0, 0, 0};
    bcf_hdr_format(vcf_header, 0, &htxt);
    std::cout << htxt.s;
    free(htxt.s);
    bcf_hdr_destroy(vcf_header);
    bcf_destroy(vcf_record);
    bcf_close(vcf);
}
void step2(CKMCFile& kmer_db, KMAP& ref_bf, BF& bf, const BF& context_bf)
{
    uint32 klen, mode, min_counter, pref_len, sign_len, min_c, counter;
    uint64 tot_kmers, max_c;
    kmer_db.Info(klen, mode, min_counter, pref_len,
                 sign_len, min_c, max_c,tot_kmers);
    CKmerAPI kmer_obj(klen);

    char context[opt::ref_k + 1];

    while (kmer_db.ReadNextKmer(kmer_obj, counter))
    {
        kmer_obj.to_string(context);
        std::transform(context, context + opt::ref_k, context, ::toupper);
        char kmer[opt::k + 1];
        strncpy(kmer, context + ((opt::ref_k - opt::k) / 2), opt::k);
        kmer[opt::k] = '\0';
        ref_bf.increment(kmer, counter);
        if (!context_bf.test_key(context))
        {
            bf.increment(kmer, counter);
        }
    }
}
void step3(std::vector<Variant> &vs, BF &bf, KMAP &ref_bf,
           const std::unordered_map<std::string, std::string> &refs)
{
    std::string last_seq_name = "";
    VB vb(opt::k, opt::error_rate);
    for (auto v = vs.begin(); v != vs.end(); ++v)
    {
        if (last_seq_name.size() == 0)
            last_seq_name = v->seq_name;
        if (!v->has_alts)
            continue;
        if (vb.empty())
        {
            vb.add_variant(*v);
            continue;
        }
        if (!vb.is_near_to_last(*v) || last_seq_name != v->seq_name)
        {
            VK_GROUP kmers = vb.extract_kmers(refs.at(last_seq_name));
            set_coverages(bf, ref_bf, vb, kmers);
            vb.genotype(opt::max_coverage);
            vb.output_variants(opt::verbose);
            vb.clear();
            if (last_seq_name != v->seq_name)
                last_seq_name = v->seq_name;
        }
        vb.add_variant(*v);
    }
    if (!vb.empty())
    {
        VK_GROUP kmers = vb.extract_kmers(refs.at(last_seq_name));
        set_coverages(bf, ref_bf, vb, kmers);
        vb.genotype(opt::max_coverage);
        vb.output_variants(opt::verbose);
        vb.clear();
    }
    std::cout.flush();
}

// ---------------------------------------------------------------------------
int main(int argc, char *argv[])
{
    auto start_t = omp_get_wtime();
    hts_set_log_level(HTS_LOG_OFF);

    parse_arguments(argc, argv);
    //omp_set_num_threads(2);

    CKMCFile kmer_db;
    if (!kmer_db.OpenForListing(opt::kmc_sample_path))
    {
        std::cerr << "ERROR: cannot open " << opt::kmc_sample_path << std::endl;
        return 1;
    }
    std::unordered_map<std::string, std::string> refs;
    BF bf;
    KMAP ref_bf;
    BF context_bf;
    std::vector<Variant> vs;
    std::vector<std::string> used_seq_names;
    std::future<void> bf_rank;
    double t0, t1;
    t0 = omp_get_wtime();
    //omp_set_nested(1);
    #pragma omp parallel num_threads(6) default(shared) //proc_bind(spread)
    {
        #pragma omp single nowait
        {
            #pragma omp task depend(out: refs)
            {
                double s = omp_get_wtime();
                refs = read_references();
                double e = omp_get_wtime();
                std::cerr<<"references read in "<<e-s<<std::endl;
            }
            #pragma omp task depend(in: refs)
            {
                double s = omp_get_wtime();
                omp_set_nested(1);
                for(auto it = refs.begin(); it != refs.end(); ++it)
                {
                    std::string id = it->first;
                    upper_sequence(refs[id]);
                }
                double e = omp_get_wtime();
                std::cerr<<"upped all sequences in:  "<<e-s<<std::endl;


            }
            #pragma omp task
            {
                double s = omp_get_wtime();
                omp_set_nested(1);
                bf = BF(opt::bf_size);
                double e = omp_get_wtime();
                std::cerr<<"bf init done in "<<e-s<<std::endl;
            }
        #pragma omp task
            {
                double s = omp_get_wtime();
                omp_set_nested(1);
                context_bf = BF(opt::bf_size);
                double e = omp_get_wtime();
                std::cerr<<"context_bf init in "<<e-s<<std::endl;
            }
        #pragma omp task
        {
                double s = omp_get_wtime();
            vs = read_variants();
            for (auto v : vs)
            {
                if(used_seq_names.empty())
                    used_seq_names.push_back(v.seq_name);
                else
                    if(used_seq_names.back() != v.seq_name)
                        used_seq_names.push_back(v.seq_name);
                
            }
                double e = omp_get_wtime();
                std::cerr<<"variants read in "<<e-s<<std::endl;
        }
        }
    }
    t1 = omp_get_wtime();
    std::cerr<<"Preprocessing done in: "<<t1-t0<<"s\n";

    compute_main_filters(vs, bf, ref_bf, refs);

    t0 = omp_get_wtime();
    
    std::cerr<<"Bf done in: "<<t0-t1<<"s\n";
    pelapsed("BF creation complete", start_t);
    
    std::cerr<<"num device: "<<omp_get_num_devices()<<std::endl;
    //bf_rank = std::async(std::launch::async, [&]() { bf.switch_mode(); });   

    //t0 = omp_get_wtime();
                double s = omp_get_wtime();
                omp_set_nested(1);
                bf.switch_mode();
                double e = omp_get_wtime();
                std::cerr<<"switch done in "<<e-s<<std::endl;
                 s = omp_get_wtime();
                compute_context_filter(bf, context_bf, refs, used_seq_names);
                context_bf.read_mode();
                 e = omp_get_wtime();
                std::cerr<<"context done in "<<e-s<<std::endl;


                /* #pragma omp parallel
    {
        #pragma omp single nowait
        {
            #pragma omp task
            {
                double s = omp_get_wtime();
                omp_set_nested(1);
                bf.switch_mode();
                double e = omp_get_wtime();
                std::cerr<<"switch done in "<<e-s<<std::endl;
            }
            #pragma omp task
            {
                double s = omp_get_wtime();
                compute_context_filter(bf, context_bf, refs, used_seq_names);
                context_bf.read_mode();
                double e = omp_get_wtime();
                std::cerr<<"context done in "<<e-s<<std::endl;
            }
                */
    
            //t1 = omp_get_wtime();
            //std::cerr << "Context_bf done in: " << t1 - t0 <<"s\n";
            
            /*std::ofstream file;
              file.open("strings2.txt");
              for(auto s: cv)
              file << s<<std::endl;
              file.close();
            */
            // context_bf.switch_mode();
            //t0 = omp_get_wtime();
            //std::cerr<<"Waiting...\n";
            //bf_rank.wait();
            //t1 = omp_get_wtime();
            //std::cerr<<"Waited for "<<t1 - t0 <<"s\n";
                //}
                //}
            
    pelapsed("Reference BF creation complete", start_t);
            
    // STEP 2: test variants present in read sample
    t0 = omp_get_wtime();
    step2(kmer_db, ref_bf, bf, context_bf);
    t1 = omp_get_wtime();
    std::cerr << "Weights done in " << t1 - t0 <<"s\n";
    pelapsed("BF weights created", start_t);

    // STEP 3: check if variants in vcf are covered enough
    print_cleaned_header();
    step3(vs, bf, ref_bf, refs);

    t0 = omp_get_wtime();
    std::cerr << "Print to cout done in " << t0 - t1 <<"s\n";

    pelapsed("Execution completed", start_t);

    return 0;
}

/*void populate_context_bf(const BF &bf, BF &context_bf,
                         const std::unordered_map<std::string, std::string> &refs,
                         const std::vector<std::string> &used_seq_names) {
    if(used_seq_names.size() == 0)
        return;
    for (long unsigned int seq = 0; seq < used_seq_names.size(); ++seq)
    {
        long unsigned int size = refs.at(used_seq_names[seq]).size();
        const char* ref = refs.at(used_seq_names[seq]).c_str();
        int pos = (opt::ref_k - opt::k) / 2;
        //#pragma omp parallel for proc_bind(spread)#pragma omp end declare target
        //std::cerr<<"num devices: "<<omp_get_num_devices()<<std::endl;
        uint64_t* cv = new uint64_t[size];
        #pragma acc parallel loop independent gang vector copyin(ref[:size], pos) copyout(cv[:size])
        for (long unsigned int p = opt::ref_k; p < size; ++p)
        {
            auto it0 = p - opt::ref_k;
            auto it1 = p - opt::ref_k + pos;
            char t[64];
            char ckmer[64];
            #pragma acc loop independent vector(64)
            for(int i =0 ; i < opt::k ; ++i)
                t[i] = ref[it1 + i];
            //cv[it0]=calc_hash(t, opt::k);
            //reverse_cmpl(ckmer, kmer, k);
            #pragma acc loop independent vector(64)
            for(int i = 0; i < opt::k; ++i) {
                ckmer[i] = RCN[(int)t[opt::k - 1 - i]];
            }
            int i = 0;
            #pragma acc loop seq private(i)
            while(i < opt::k && t[i] == ckmer[i])
                ++i;
            uint64_t hash[2];
            if(i == opt::k || t[i] < ckmer[i])
                MurmurHash3_x64_128(t, (int)opt::k, 0, reinterpret_cast<void *>(&hash));
            else
                MurmurHash3_x64_128(ckmer, (int)opt::k, 0, reinterpret_cast<void *>(&hash));
            cv[it0] = hash[0];	
            
            //if (bf.test_key(std::string_view(&reference[it1], &reference[it12])))
            
            //if(bf.test_key(std::string_view(&reference[it1], opt::k)))
            //  context_bf.add_key(std::string_view(&reference[it0], opt::ref_k));
        }
        for(long unsigned int i = 0; i < size; ++i)
            if(bf.ttest(cv[i]))
                
            delete[] cv;
    }
}

*/
