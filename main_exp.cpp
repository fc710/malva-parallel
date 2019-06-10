
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
#include <chrono>
#include <fstream>
#include <iostream>
#include <utility>

#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include <math.h>
#include <zlib.h>

//#include "tbb/tbb.h"
#include "tbb/concurrent_queue.h"
#include "tbb/concurrent_unordered_map.h"
#include "tbb/concurrent_vector.h"
#include "tbb/parallel_for.h"
#include "tbb/pipeline.h"
#include "tbb/spin_mutex.h"
#include <omp.h>

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

auto start_t = std::chrono::high_resolution_clock::now();

void pelapsed(const std::string &s = "") {
		auto now_t = std::chrono::high_resolution_clock::now();
		std::cerr << "[malva-geno/" << s << "] Time elapsed "
				  << std::chrono::duration_cast<std::chrono::milliseconds>(
						  now_t -start_t).count()/1000 << "s" << std::endl;
}

KSEQ_INIT(gzFile, gzread)

std::unordered_map<std::string, std::string> read_references()
{
		gzFile fasta_in = gzopen(opt::fasta_path.c_str(), "r");
		kseq_t *reference = kseq_init(fasta_in);
		// int l;
		std::unordered_map<std::string, std::string> refs;
		while (kseq_read(reference) >= 0)
		{
				std::string id = reference->name.s;
				if (id.compare(0, 3, "chr") == 0) {
						id = id.substr(3);
				}
				std::string seq(reference->seq.s);
				std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
				refs[id] = seq;
		}

		kseq_destroy(reference);
		gzclose(fasta_in);
		return refs;
}

std::vector<Variant> read_variants()
{
		htsFile *vcf = bcf_open(opt::vcf_path.c_str(), "r");
		bcf_hdr_t *vcf_header = bcf_hdr_read(vcf);
		bcf1_t *vcf_record = bcf_init();
		std::vector<Variant> vs;
		while (bcf_read(vcf, vcf_header, vcf_record) == 0)
		{
				bcf_unpack(vcf_record, BCF_UN_STR);
				vs.emplace_back(vcf_header, vcf_record, opt::pop);
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
		long unsigned int size = kmers.size();
		for (long unsigned int i = 0; i < size; ++i)
		{
				// For each variant
				const auto &v = kmers.at(i);
				for (const auto &p : v)
						// For each allele of the variant/
						for (const auto &Ks : p.second)
								// For each list of kmers of the allele
								for (const auto &kmer : Ks)
										// For each kmer in the kmer list
										if (p.first == 0)
												ref_bf.add_key(kmer);
										else
												bf.add_key(kmer);
		}
}

void populate_bf_refbf(std::vector<Variant> &vs, BF &bf, KMAP &ref_bf,
					   const std::unordered_map<std::string, std::string> &refs)
{
#pragma omp parallel
		{
#pragma omp single
				{
						std::string last_seq_name = "";
						VB vb(opt::k, opt::error_rate);
						for (decltype(vs.size()) i = 0; i < vs.size(); ++i)
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
#pragma omp task firstprivate(vb, last_seq_name)
										{
												VK_GROUP kmers = vb.extract_kmers(refs.at(last_seq_name));
												add_kmers_to_bf(bf, ref_bf, kmers);
										}
										vb.clear();
										if (last_seq_name != v.seq_name)
												last_seq_name = v.seq_name;
								}
								vb.add_variant(v);
						}
						if (!vb.empty()) {
								VK_GROUP kmers = vb.extract_kmers(refs.at(last_seq_name));
								add_kmers_to_bf(bf, ref_bf, kmers);
						}
				}
		}
}

void populate_context_bf(const BF &bf, BF &context_bf,
						 const std::unordered_map<std::string, std::string> &refs,
						 std::vector<std::string> &used_seq_names) {
#pragma omp parallel
		{
				for (long unsigned int seq = 0; seq < used_seq_names.size(); ++seq)
				{
						std::string reference = refs.at(used_seq_names[seq]);
						int pos = (opt::ref_k - opt::k) / 2;
						auto size = reference.size();

#pragma omp for // collapse(2)
						for (long unsigned int p = opt::ref_k; p < size; ++p)
						{
								auto it0 = p - opt::ref_k;
								auto it1 = p - opt::ref_k + pos;
								// auto it02 = p;
								auto it12 = p - pos;
								if (bf.test_key(std::string(&reference[it1], &reference[it12])))
								{
										context_bf.add_key(std::string(&reference[it0], &reference[p]));
								}
						}
				}
		}
}

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
		hts_set_log_level(HTS_LOG_OFF);

		parse_arguments(argc, argv);

		// STEP 0: open and check input files
		CKMCFile kmer_db;
		if (!kmer_db.OpenForListing(opt::kmc_sample_path))
		{
				std::cerr << "ERROR: cannot open " << opt::kmc_sample_path << std::endl;
				return 1;
		}
		// References are stored in a map
		std::unordered_map<std::string, std::string> refs;
		BF bf;
		KMAP ref_bf;
		BF context_bf;

		std::vector<Variant> vs;
		// tbb::concurrent_vector<Variant> vs;
		std::vector<std::string> used_seq_names;
		std::future<void> bf_rank;
		double t0, t1;
		// STEP 1: add VCF kmers to bloom filter
#pragma omp parallel
		{
#pragma omp single
				{
#pragma omp task
						{
								refs = read_references();
						}
#pragma omp task
						{
								bf = BF(opt::bf_size);
						}
#pragma omp task
						{
								context_bf = BF(opt::bf_size);
						}
#pragma omp task
						{
								vs = read_variants();
								for (const auto &v : vs)
										used_seq_names.push_back(v.seq_name);
								sort(used_seq_names.begin(), used_seq_names.end());
								used_seq_names.erase(unique(used_seq_names.begin(), used_seq_names.end()),
													 used_seq_names.end());
						}
				}
		}
		populate_bf_refbf(vs, bf, ref_bf, refs);
		bf_rank = std::async(std::launch::async, [&]() { bf.switch_mode(); });
		pelapsed("BF creation complete");
		t0 = omp_get_wtime();
		populate_context_bf(bf, context_bf, refs, used_seq_names);
		// std::transform(reference.begin(), reference.end(), reference.begin(),
		// ::toupper); 	std::string ref_ksub(reference, (opt::ref_k - opt::k) / 2,
		//opt::k); 	std::string context(reference, 0, opt::ref_k);
		t1 = omp_get_wtime();
		std::cerr << "ctx creation took: " << t1 - t0 << std::endl;

		/*std::ofstream file;
		  file.open("strings2.txt");
		  for(auto s: cv)
		  file << s<<std::endl;
		  file.close();
		*/
		pelapsed("Reference BF creation complete");
		// context_bf.switch_mode();
		context_bf.read_mode();

		bf_rank.wait();
		// STEP 2: test variants present in read sample
		double t8 = omp_get_wtime();
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
		double t9 = omp_get_wtime();
		std::cerr << "weights creation in " << t9 - t8 << std::endl;
		pelapsed("BF weights created");

		// STEP 3: check if variants in vcf are covered enough
		print_cleaned_header();
		step3(vs, bf, ref_bf, refs);

		double t10 = omp_get_wtime();
		std::cerr << "last step took " << t10 - t9 << std::endl;

		pelapsed("Execution completed");

		return 0;
}
/*		tbb::parallel_for(tbb::blocked_range<int>(opt::ref_k, size),
		[&](const tbb::blocked_range<int>& r){
		int b =r.begin();
		int e =r.end();
		for (int p = b,it0 = p - opt::ref_k,it1 = p - opt::ref_k + pos,
		it02 = p ,it12 = p - pos; p != e; ++p,++it0,++it1,++it02,++it12)
		{

		if(bf.test_key(std::string(&reference[it1], &reference[it12]))){
		context_bf.add_key(std::string(&reference[it0],
		&reference[it02]));
		}
		}
		}
		);
*/
