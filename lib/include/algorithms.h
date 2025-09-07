#ifndef ALGORITHMS_H
#define ALGORITHMS_H


#ifdef USE_OPENMP
#include <omp.h>
#endif

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <omp.h>
#include <cstdio>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <unordered_map>


#ifdef __cplusplus
extern "C" {
#endif

extern void hash3_vector(std::vector<uint32_t>& freq ,const std::string& input_file, const uint32_t len_, const bool perf_collect);
extern void naive_vector(std::vector<uint32_t>& freq ,const std::string& input_file, const uint32_t len_, const bool perf_collect);
extern void base_naive(std::vector<uint32_t>& freq ,const std::string& input_file, const uint32_t len_);
extern void rabin_karp_rolling_hash_vector(std::vector<uint32_t>& freq, const std::string& input_file, const uint32_t len_, const bool perf_collect);
extern void rabin_karp_SWAR_vector(std::vector<uint32_t>& freq, const std::string& input_file, const uint32_t len_, const bool perf_collect);
extern void hash3_scalar(std::vector<uint32_t>& freq ,const std::string& input_file, const uint32_t len_, const bool perf_collect);
extern void naive_scalar(std::vector<uint32_t>& freq ,const std::string& input_file, const uint32_t len_, const bool perf_collect);
extern void rabin_karp_rolling_hash_scalar(std::vector<uint32_t>& freq ,const std::string& input_file, const uint32_t len_, const bool perf_collect);
extern void rabin_karp_SWAR_scalar(std::vector<uint32_t>& freq ,const std::string& input_file, const uint32_t len_, const bool perf_collect) ;

#ifdef __cplusplus
}
#endif

#endif // ALGORITHMS_H