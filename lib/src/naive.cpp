#include "algorithms.h"
#include "common_defs.h"
#include <Kokkos_Core.hpp>
#include <Kokkos_SIMD.hpp>

typedef Kokkos::Experimental::native_simd<int8_t> vs8;

extern void naive(std::vector<uint32_t>& freq, 
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
    size_t size = data_.size();
    int len = len_;

    // Determine SIMD vector length
    constexpr size_t VECTOR_LENGTH = vs8::size();
    uint8_t cycles = len / VECTOR_LENGTH;
    uint8_t leftover = len % VECTOR_LENGTH;

    // Prepare symbol mapping
    std::unordered_map<int8_t, int8_t> symbols_code {
        {'A', int8_t(0)}, {'C', int8_t(1)}, {'G', int8_t(2)}, {'T', int8_t(3)}
    };

    // Create Kokkos views
    Kokkos::View<int8_t*> data("data", size);
    for(int i = 0; i < size; ++i) {
        data(i) = symbols_code[data_[i]];
    }

    Kokkos::View<uint32_t*> freq_view("freq_view", size);

    timer.reset();

    // Parallel execution with Kokkos
    Kokkos::parallel_for("naive_string_match", size - len + 1, KOKKOS_LAMBDA(const int i) {
        uint32_t res = 0;
        
        for (int j = 0; j < size - len + 1; ++j) {
            int is_eq = 1;
            
            // Vectorized comparison
            for (int k = 0; k < cycles && is_eq; ++k) {
                vs8 vpattern_1(&data(i + k * VECTOR_LENGTH), Kokkos::Experimental::element_aligned_tag());
                vs8 vpattern_2(&data(j + k * VECTOR_LENGTH), Kokkos::Experimental::element_aligned_tag());
                
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
            }
            
            res += is_eq;
        }
        
        freq_view[i] = res;
    });

    Kokkos::fence();
    double stop = timer.seconds();

    // Copy results back to host
    freq.resize(size);
    auto freq_host = Kokkos::create_mirror_view(freq_view);
    Kokkos::deep_copy(freq_host, freq_view);
    
    for (size_t i = 0; i < size; ++i) {
        freq[i] = freq_host[i];
    }

    if (perf_collect) {
        std::cout << freq.size() << " " << len << " " << stop << "\n";
    }
}
    Kokkos::finalize();
    return;
}