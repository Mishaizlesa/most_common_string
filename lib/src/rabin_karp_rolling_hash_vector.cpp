#include <sycl/sycl.hpp>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <unordered_map>
#include "common_defs.h"

extern "C" void rabin_karp_rolling_hash_vector(std::vector<uint32_t>& freq, const std::string& input_file, 
                                       const uint32_t len_, const bool perf_collect) {
    std::ifstream fin(input_file);
    std::string data_str;
    fin >> data_str;
    uint64_t size = data_str.size();
    int len = len_;

    freq.resize(size, 0);

    std::vector<uint8_t> data(size);
    std::unordered_map<char, uint8_t> mapSymbToCode = {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};
    
    for (uint64_t i = 0; i < size; ++i) {
        data[i] = mapSymbToCode[data_str[i]];
    }

    auto start = std::chrono::high_resolution_clock::now();

        sycl::queue q(sycl::cpu_selector_v);
        
        const size_t cycles = len / VEC_LEN;
        const size_t leftover = len % VEC_LEN;
        
        sycl::buffer<uint8_t> data_buf(data.data(), size);
        sycl::buffer<uint32_t> freq_buf(freq.data(), size);

        q.submit([&](sycl::handler& h) {
            auto data_acc = data_buf.get_access<sycl::access::mode::read>(h);
            auto freq_acc = freq_buf.get_access<sycl::access::mode::write>(h);

            h.parallel_for(sycl::range<1>(size - len + 1), [=](sycl::id<1> idx) {
                uint64_t i = idx[0];
                uint32_t res = 0;
                
                uint64_t hash_pattern = 0;
                uint64_t p = 2;
                uint64_t powmod = 1;
                
                for (int j = 0; j < len; ++j) {
                    hash_pattern = hash_pattern * p + data_acc[i + j];
                    powmod = powmod * p;
                }
                
                uint64_t hash_text = 0;
                for (int j = 0; j < len; ++j) {
                    hash_text = hash_text * p + data_acc[j];
                }
                
                for (int j = 0; j < size - len + 1; ++j) {
                    if (hash_text == hash_pattern) {
                        bool is_eq = true;
                        
                        for (size_t k = 0; k < cycles && is_eq; ++k) {
                            const size_t offset1 = i + k * VEC_LEN;
                            const size_t offset2 = j + k * VEC_LEN;
                            
                            sycl::vec<uint8_t, VEC_LEN> vec_pattern, vec_text;
                            
                            vec_pattern.load(0, &data_acc[offset1]);
                            vec_text.load(0, &data_acc[offset2]);
                            
                            auto mask = (vec_pattern == vec_text);
                            
                            is_eq = sycl::all(mask);
                            if (!is_eq) break;
                        }
                        
                        // Handle leftover elements
                        for (size_t k = 0; k < leftover && is_eq; ++k) {
                            const size_t pos1 = i + cycles * VEC_LEN + k;
                            const size_t pos2 = j + cycles * VEC_LEN + k;
                            if (data_acc[pos1] != data_acc[pos2]) {
                                is_eq = false;
                                break;
                            }
                        }
                        
                        if (is_eq) {
                            res++;
                        }
                    }
                    
                    if (j < size - len) {
                        hash_text = (hash_text * p - data_acc[j] * powmod + data_acc[j + len]);
                    }
                }
                
                freq_acc[i] = res;
            });
        });
        
        q.wait();


    auto stop = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(stop - start).count();

    if (perf_collect) {
        std::cout << size << " " << len_ << " " << elapsed << "\n";
    }
}