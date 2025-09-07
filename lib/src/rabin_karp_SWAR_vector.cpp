#include <sycl/sycl.hpp>
#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include <iostream>
#include <chrono>
#include "common_defs.h"

extern "C" void rabin_karp_SWAR_vector(std::vector<uint32_t>& freq, const std::string& input_file, const uint32_t len_, const bool perf_collect) {
    std::ifstream fin(input_file);
    std::string data_;
    fin >> data_;
    uint64_t size = data_.size();
    int len = len_;
    int32_t cycles = len / VEC_LEN;
    int32_t leftover = len % VEC_LEN;

    freq.resize(size, 0);
    std::vector<uint8_t> data(size + 65, 5);
    std::unordered_map<char, int8_t> mapSymbToCode = {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};

    for (size_t i = 0; i < size; ++i) {
        data[i] = mapSymbToCode[data_[i]];
    }

    auto start = std::chrono::high_resolution_clock::now();

    sycl::queue q(sycl::default_selector_v);

    sycl::buffer<uint8_t> data_buf(data.data(), sycl::range<1>(size + 65));
    sycl::buffer<uint32_t> freq_buf(freq.data(), sycl::range<1>(size));

    q.submit([&](sycl::handler& cgh) {
        auto data_acc = data_buf.get_access<sycl::access::mode::read>(cgh);
        auto freq_acc = freq_buf.get_access<sycl::access::mode::write>(cgh);

        using vec16 = sycl::vec<uint8_t, 16>;

        cgh.parallel_for(sycl::range<1>(size - len + 1), [=](sycl::id<1> idx) {
            int i = idx[0];
            uint32_t res = 0;

            uint8_t pattern_first1 = data_acc[i];
            uint8_t pattern_first2 = data_acc[i + 1];
            uint8_t pattern_last1 = data_acc[i + len - 1];
            uint8_t pattern_last2 = data_acc[i + len - 2];

            for (uint64_t j = 0; j < size - len + 1; j += VEC_LEN_SWAR) {
                uint64_t bitmask = 0;

                for (uint64_t k = 0; k < VEC_LEN_SWAR; k += VEC_LEN) {
                    vec16 vfirst_sym1, vfirst_sym2, vlast_sym1, vlast_sym2;
                    uint8_t temp1[16], temp2[16], temp3[16], temp4[16];
                    for (int m = 0; m < 16; ++m) {
                        temp1[m] = data_acc[j + k + m];
                        temp2[m] = data_acc[j + k + 1 + m];
                        temp3[m] = data_acc[j + k + len - 1 + m];
                        temp4[m] = data_acc[j + k + len - 2 + m];
                    }
                    vfirst_sym1.load(0, temp1);
                    vfirst_sym2.load(0, temp2);
                    vlast_sym1.load(0, temp3);
                    vlast_sym2.load(0, temp4);

                    vec16 vpattern_first1(pattern_first1);
                    vec16 vpattern_first2(pattern_first2);
                    vec16 vpattern_last1(pattern_last1);
                    vec16 vpattern_last2(pattern_last2);

                    vec16 eq_first1 = vfirst_sym1 ^ vpattern_first1;
                    vec16 eq_first2 = vfirst_sym2 ^ vpattern_first2;
                    vec16 eq_last1 = vlast_sym1 ^ vpattern_last1;
                    vec16 eq_last2 = vlast_sym2 ^ vpattern_last2;
                    vec16 veq_first = eq_first1 | eq_first2;
                    vec16 veq_last = eq_last1 | eq_last2;
                    vec16 veq = veq_first | veq_last;

                    for (int m = 0; m < 16; ++m) {
                        if (veq[m] == 0) {
                            bitmask |= (1ULL << (k + m));
                        }
                    }
                }

                for (uint64_t t = 0; t < VEC_LEN_SWAR; ++t) {
                    if (bitmask & (1ULL << t)) {
                        bool is_eq = true;
                        const uint8_t* pattern_1 = &data_acc[i];
                        const uint8_t* pattern_2 = &data_acc[j + t];

                        for (int k = 0; k < cycles; ++k) {
                            vec16 vpattern_1, vpattern_2;
                            vpattern_1.load(0, &pattern_1[k * VEC_LEN]);
                            vpattern_2.load(0, &pattern_2[k * VEC_LEN]);
                            auto eq_mask = (vpattern_1 == vpattern_2);
                            is_eq = sycl::all(eq_mask);
                            if (!is_eq) break;
                        }

                        for (int k = 0; k < leftover && is_eq; ++k) {
                            if (pattern_1[cycles * VEC_LEN + k] != pattern_2[cycles * VEC_LEN + k]) {
                                is_eq = false;
                                break;
                            }
                        }

                        res += is_eq ? 1 : 0;
                    }
                }
            }
            freq_acc[i] = res;
        });
    }).wait(); 

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() / 1000000.0;

    if (perf_collect) {
        std::cout << freq.size() << " " << len << " " << duration << "\n";
    }
}