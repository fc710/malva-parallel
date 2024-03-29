/**
 * MALVA - genotyping by Mapping-free ALternate-allele detection of known VAriants
 * Copyright (C) 2019  Giulia Bernardini, Luca Denti, Marco Previtali
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

#ifndef _KMAP_HPP_
#define _KMAP_HPP_

#include <string>
#include <unordered_map>
#include <algorithm>
#include <cstring>
// static const char RCN[128] = {
//     0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   //  0
//     0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 10
//     0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 20
//     0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 30
//     0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 40
//     0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 50
//     0,   0,   0, 0,   0,   'T', 0,   'G', 0,   0,   // 60
//     0,   'C', 0, 0,   0,   0,   0,   0,   'N', 0,   // 70
//     0,   0,   0, 0,   'A', 0,   0,   0,   0,   0,   // 80
//     0,   0,   0, 0,   0,   0,   0,   'T', 0,   'G', // 90
//     0,   0,   0, 'G', 0,   0,   0,   0,   0,   0,   // 100
//     'N', 0,   0, 0,   0,   0,   'A', 0,   0,   0,   // 110
//     0,   0,   0, 0,   0,   0,   0,   0              // 120
// };
#include "bloom_filter.hpp"
struct KMAP {
  std::unordered_map<std::string, int> kmers;
  std::unordered_map<std::string, int> _times;

  KMAP() {}

	static char _compl(const char &c) { return opt::RCN[(int)c]; }

  std::string canonical(const char* kmer) {
    uint k = strlen(kmer);
    char ckmer[k + 1];
    strcpy(ckmer, kmer);
    std::transform(ckmer, ckmer + k, ckmer, _compl);
    std::reverse(ckmer, ckmer + k);
    if (strcmp(kmer, ckmer) < 0)
      memmove(ckmer, kmer, k);
    std::string kmer_string (ckmer);
    return kmer_string;
  }
	std::string canonical(const std::string& kmer) const {
		std::string ckmer(kmer);
		std::transform(ckmer.begin(), ckmer.end(), ckmer.begin(), _compl);
		std::reverse(ckmer.begin(), ckmer.end());
		if(kmer.compare(ckmer) < 0)
			return kmer;
		else
			return ckmer;	
	}
	std::string _reverse_cmpl(std::string_view kmer) const
	{
		std::string ckmer(kmer);
		int size = ckmer.size();
        #pragma ivdep
		for(int i = 0; i < size; ++i)
				ckmer[i] = opt::RCN[(int)ckmer[i]];
		std::reverse(ckmer.begin(), ckmer.end());
			return ckmer;
		
	}



  bool test_key(const char* kmer) {
    std::string ckmer = canonical(kmer);
    if(kmers.find(ckmer) == kmers.end())
      return false;
    else
      return true;
  }

  void add_key(const char* kmer) {
    std::string ckmer = canonical(kmer);
    kmers[ckmer] = 0;
  }
	void add_key(std::string_view kmer){
		std::string ckmer = _reverse_cmpl(kmer);
		if(kmer.compare(ckmer) > 0)
			kmer = ckmer;
		#pragma omp critical 
		kmers[(std::string)kmer] = 0;		
	}

	void add_key(std::string kmer) {
		std::string ckmer = canonical(kmer);
		kmers[ckmer] = 0;
	}


  void increment(const char* kmer, int counter) {
    std::string ckmer = canonical(kmer);
    if(kmers.find(ckmer) != kmers.end()) {
      uint32_t new_value = kmers[ckmer] + counter;
      kmers[ckmer] = new_value < 250 ? new_value : 250;
      ++_times[ckmer];
    }
  }

  void increment_with_average(const char* kmer, int counter) {
    std::string ckmer = canonical(kmer);
    if(kmers.find(ckmer) != kmers.end()) {
      uint32_t new_value = (kmers[ckmer] * _times[ckmer] + counter) / (_times[ckmer]+1);
      kmers[ckmer] = new_value < 250 ? new_value : 250;
      ++_times[ckmer];
    }
  }

  int get_count(const char* kmer) {
    std::string ckmer = canonical(kmer);
    if(kmers.find(ckmer) != kmers.end())
      return kmers[ckmer];
    else
      return 0;
  }

  int get_times(const char* kmer) {
    std::string ckmer = canonical(kmer);
    if(kmers.find(ckmer) != kmers.end())
      return _times[ckmer];
    else
      return 0;
  }
};

#endif
