#include "algorithms.h"
#include "common_defs.h"
#include <Kokkos_Core.hpp>
#include <Kokkos_SIMD.hpp>

typedef Kokkos::Experimental::native_simd<uint64_t> vu64;
typedef Kokkos::Experimental::native_simd<uint8_t> vu8;

extern void rabin_karp_rolling_hash_scalar(std::vector<uint32_t>& freq, 
                                  const std::string& input_file, 
                                  const uint32_t len_, 
                                  const bool perf_collect) {
    // Initialize Kokkos
    Kokkos::initialize();
    {
        Kokkos::Timer timer;
    std::ifstream fin(input_file);
    std::string data_;
    fin >> data_;
    uint64_t size = data_.size();
    int len = len_;

    constexpr size_t VECTOR_LENGTH = vu8::size();
    int32_t cycles = len / VECTOR_LENGTH;
    int32_t leftover = len % VECTOR_LENGTH;

    freq.resize(size, 0);
    
    // Prepare data
    Kokkos::View<uint8_t*> data("data", size + len);
    std::unordered_map<int8_t, int8_t> mapSymbToCode = {
        {'A', (int8_t)0}, {'C', (int8_t)1}, {'G', (int8_t)2}, {'T', (int8_t)3}
    };

    // Initialize data on host
    auto data_host = Kokkos::create_mirror_view(data);
    for (int i = 0; i < size; ++i) {
        data_host[i] = mapSymbToCode[data_[i]];
    }
    for (int i = size; i < size + len; ++i) {
        data_host[i] = 7;
    }
    Kokkos::deep_copy(data, data_host);

    // Create Kokkos view for frequencies
    Kokkos::View<uint32_t*> freq_view("freq_view", size);
    
    timer.reset();

    // Parallel execution with Kokkos
    Kokkos::parallel_for("rabin_karp", size - len + 1, KOKKOS_LAMBDA(const int i) {
        int res = 0;
        uint64_t hash_pattern = 0;
        uint64_t hash_text = 0;
        uint64_t p = 2;
        uint64_t powmod = 1;
        
        // Compute initial hashes
        for (int j = 0; j < len; ++j) {
            hash_pattern = hash_pattern * p + data[i + j];
            hash_text = hash_text * p + data[j];
            powmod = powmod * p;
        }
        
        // Search for matches
        for (int j = 0; j < size - len + 1; ++j) {
            if (hash_text == hash_pattern) {
                int is_eq = 1;


                for(int k = 0; k < len; ++k)
                {
                    if (data[j + k] != data[i + k])
                    {
                        is_eq = 0;
                        break;
                    }
                }
                
                // Vectorized comparison
                /*for (int k = 0; k < cycles && is_eq; ++k) {
                    vu8 vpattern_1(&data[i + k * VECTOR_LENGTH], Kokkos::Experimental::element_aligned_tag{});
                    vu8 vpattern_2(&data[j + k * VECTOR_LENGTH], Kokkos::Experimental::element_aligned_tag{});
                    
                    is_eq = Kokkos::Experimental::all_of(vpattern_1 == vpattern_2);
                    if (is_eq == 0) {
                        break;
                    }
                }
                
                // Leftover elements
                for (int k = 0; k < leftover && is_eq; ++k) {
                    if (data[i + cycles * VECTOR_LENGTH + k] != 
                        data[j + cycles * VECTOR_LENGTH + k]) {
                        is_eq = 0;
                        break;
                    }
                }*/
                
                res += is_eq;
            }
            
            // Update rolling hash
            if (j < size - len) {
                hash_text = (hash_text * p - data[j] * powmod + data[j + len]);
            }
        }
        
        freq_view[i] = res;
    });
    
    Kokkos::fence();
    double stop = timer.seconds();
    
    // Copy results back to host
    auto freq_host = Kokkos::create_mirror_view(freq_view);
    Kokkos::deep_copy(freq_host, freq_view);
    
    // Copy to output vector
    for (int i = 0; i < size; ++i) {
        freq[i] = freq_host[i];
    }
    
    if (perf_collect) {
        std::cout << freq.size() << " " << len << " " << stop << "\n";
    }
}
    Kokkos::finalize();
    return;
}