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

#ifndef _BLOOM_FILTER_HPP
#define _BLOOM_FILTER_HPP

#include <algorithm>
#include <array>
#include <cstring>
#include "sdsl/bit_vectors.hpp"
#include <string_view>
//#include "tbb/spin_mutex.h"

#include "MurmurHash3.hpp"
#include "kmc_api/kmc_file.h"

// using namespace std;
using namespace sdsl;

namespace opt{
	const char RCN[128] = {
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   //  0
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 10
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 20
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 30
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 40
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 50
    0,   0,   0, 0,   0,   'T', 0,   'G', 0,   0,   // 60
    0,   'C', 0, 0,   0,   0,   0,   0,   'N', 0,   // 70
    0,   0,   0, 0,   'A', 0,   0,   0,   0,   0,   // 80
    0,   0,   0, 0,   0,   0,   0,   'T', 0,   'G', // 90
    0,   0,   0, 'G', 0,   0,   0,   0,   0,   0,   // 100
    'N', 0,   0, 0,   0,   0,   'A', 0,   0,   0,   // 110
    0,   0,   0, 0,   0,   0,   0,   0              // 120
	};
}

class BF {

private:
#pragma omp declare simd
		static inline char _compl(char c) { return opt::RCN[(int)c]; }

  void _canonical(const char *kmer, char *ckmer, const int &k) const {
	strcpy(ckmer, kmer);
    std::transform( ckmer, ckmer + k, ckmer, _compl);
    std::reverse(ckmer, ckmer + k);
    if (strcmp(kmer, ckmer) < 0)
      memmove(ckmer, kmer, k);
  }

	std::string _reverse_cmpl(std::string_view kmer) const {
		std::string ckmer(kmer);
		int size = ckmer.size();
        #pragma ivdep
		for(int i = 0; i < size; ++i)
				ckmer[i] = opt::RCN[(int)ckmer[i]];
		std::reverse(ckmer.begin(), ckmer.end());
			return ckmer;
		
	}
	uint64_t _get_hash(const char *kmer) const {
		uint k = strlen(kmer);
		char ckmer[k + 1];
		_canonical(kmer, ckmer, k);
		std::array<uint64_t, 2> hashes;
		MurmurHash3_x64_128(ckmer, k, 0, reinterpret_cast<void *>(&hashes));
		return hashes[0];
	}
	
	uint64_t _get_hash(std::string_view kmer) const {
		std::string ckmer = _reverse_cmpl(kmer);
		std::array<uint64_t, 2> hashes;
		if(kmer.compare(ckmer) < 0)
				MurmurHash3_x64_128(&kmer.at(0), (uint)kmer.size(), 0, reinterpret_cast<void *>(&hashes));
		else
				MurmurHash3_x64_128(&ckmer.at(0), (uint)ckmer.size(), 0, reinterpret_cast<void *>(&hashes));
		return hashes[0];		
   }
public:
	friend void swap(BF& first, BF& second) {
		using std::swap;
		swap(first._mode, second._mode);
		swap(first._size, second._size);
		swap(first._bf, second._bf);
		swap(first._brank, second._brank);
		swap(first._counts, second._counts);
		swap(first._times, second._times);
	}
	BF()  = default;
	explicit BF(const size_t size) : _mode(false), _bf(size, 0) { _size = size; }
	BF(const BF& other) : _mode(other._mode), _size(other._size), _bf(other._bf), _brank(other._brank),
						_counts(other._counts), _times(other._times) {
	  std::cerr<<"bf copy constr"<<std::endl;
	}
	/*BF(BF&& other) :_mode(std::move(other._mode)), _size(std::move(other._size)),
					_bf(std::move(other._bf)), _brank(std::move(other._brank)),
					_counts(std::move(other._counts)), _times(std::move(other._times)) {
	  std::cout<<"bf move constr"<<std::endl;
	}
	*/
	BF(BF&& other) :BF() {
		swap(*this, other);
		std::cerr<<"bf move constr"<<std::endl;
	}

	BF& operator=(BF other) {
		swap(*this, other);
		return *this;
	}
	
	~BF() = default;

  void add_key(const char *kmer)  {
    uint64_t hash = _get_hash(kmer);
    _bf[hash % _size] = 1;
  }
	//#pragma omp declare simd
	void add_key(std::string_view kmer) {
		uint64_t hash = _get_hash(kmer);
		#pragma omp critical
		_bf[hash % _size] = 1;
	}
/*	void add_key(const std::string& kmer, tbb::spin_mutex& mtx){
		uint64_t hash = _get_hash(kmer);
		tbb::spin_mutex::scoped_lock lock(mtx);
		_bf[hash % _size] = 1;
	}
*/
  void add_refkey(const char *kmer) {
    uint64_t hash = _get_hash(kmer);
    _bf[hash % _size] = 0;
  }

  bool test_key(const char *kmer) const {
    uint64_t hash = _get_hash(kmer);
    return _bf[hash % _size];
  }
	//#pragma omp declare simd
	bool test_key(std::string_view kmer)  const {
		uint64_t hash = _get_hash(kmer);
		return _bf[hash % _size];
	}

  int get_times(const char *kmer) const {
    uint64_t hash = _get_hash(kmer);
    size_t bf_idx = hash % _size;
    size_t cnts_idx = _brank(bf_idx);
    return _times[cnts_idx];
  }

  void switch_mode() {
    _mode = true;
    _brank = rank_support_v<1>(&_bf);
    _counts = int_vector<8>(_brank(_size), 0, 8);
    _times = int_vector<8>(_brank(_size), 0, 8);
  }
	void read_mode() {
		_mode = true;
	}

  bool increment(const char *kmer, const uint32 counter) {
    if (!_mode)
      return false;
    uint64_t hash = _get_hash(kmer);
    size_t bf_idx = hash % _size;
    if (_bf[bf_idx]) {
      size_t cnts_idx = _brank(bf_idx);
      uint32 new_value = _counts[cnts_idx] + counter;
      _counts[cnts_idx] = new_value < 250 ? new_value : 250;
      ++_times[cnts_idx];
    }
    return true;
  }

    bool increment_with_average(const char *kmer, const uint32 counter) {
    if (!_mode)
      return false;
    uint64_t hash = _get_hash(kmer);
    size_t bf_idx = hash % _size;
    if (_bf[bf_idx]) {
      size_t cnts_idx = _brank(bf_idx);
      uint32 new_value = (_counts[cnts_idx] * _times[cnts_idx] + counter) / (_times[cnts_idx]+1);
      _counts[cnts_idx] = new_value < 250 ? new_value : 250;
      ++_times[cnts_idx];
    }
    return true;
  }

  uint8_t get_count(const char *kmer) const {
    if (_mode) {
      uint64_t hash = _get_hash(kmer);
      size_t bf_idx = hash % _size;
      if (_bf[bf_idx])
        return _counts[_brank(bf_idx)];
    }
    return 0;
  }

private:
  // const BF &operator=(const BF &other) { return *this; }
  //const BF &operator=(const BF &&other) { return *this; }
	// BF &operator=(const BF &other) {return *this;}
	//BF &operator=(const BF &&other) {return *this;}

  bool _mode; // false = write, true = read
  size_t _size;
  bit_vector _bf;
  rank_support_v<1> _brank;
  int_vector<8> _counts;
  int_vector<8> _times;
};

#endif
