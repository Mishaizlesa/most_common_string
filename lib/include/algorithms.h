#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include "common_defs.h"

extern void hash3(std::vector<uint32_t>& freq ,const std::string& input_file, const uint32_t len_, const bool perf_collect);
extern void naive(std::vector<uint32_t>& freq ,const std::string& input_file, const uint32_t len_, const bool perf_collect);
extern void base_naive(std::vector<uint32_t>& freq ,const std::string& input_file, const uint32_t len_);
extern void rabin_karp(std::vector<uint32_t>& freq ,const std::string& input_file, const uint32_t len_, const bool perf_collect);

#endif // ALGORITHMS_H