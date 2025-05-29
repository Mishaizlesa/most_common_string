#include "algorithms.h"
#include "common_defs.h"
#include <Kokkos_Core.hpp>
#include <Kokkos_SIMD.hpp>

typedef long long ll;

extern void hash3(std::vector<uint32_t>& freq, const std::string& input_file, 
                 const uint32_t len_, const bool perf_collect) {
    // Initialize Kokkos
    Kokkos::initialize();
    {
        Kokkos::Timer timer;
    std::ifstream fin(input_file);
    std::string data_;
    fin >> data_;
    ll size = data_.size();
    int len = len_;

    std::unordered_map<int8_t, int8_t> symbols_code {
        {'A', uint8_t(0)}, {'C', uint8_t(1)}, {'G', uint8_t(2)}, {'T', uint8_t(3)}
    };

    // Convert data to Kokkos view
    Kokkos::View<int8_t*> data("data", size);
    for(int i = 0; i < size; ++i) {
        data(i) = symbols_code[data_[i]];
    }

    freq.resize(size);
    Kokkos::View<uint32_t*> freq_view("freq_view", size);

    // Create atomic counters for performance metrics
    Kokkos::View<uint64_t*> num_of_coll("num_of_coll", 1);
    Kokkos::View<uint64_t*> num_of_comp("num_of_comp", 1);

    // Main parallel loop
    timer.reset();
    Kokkos::parallel_for("hash3_main", Kokkos::RangePolicy<>(0, size - len + 1),
        KOKKOS_LAMBDA(const int i) {
            int res = 0;
            int sh1;
            int32_t shift[64];
            
            // Initialize shift values
            for(int k = 0; k < 64; ++k) {
                shift[k] = len - 2;
            }

            // Precompute shift values
            for(int j = 2; j < len - 1; ++j) {
                int ind = data(i + j - 2) * 16 + data(i + j - 1) * 4 + data(i + j);
                shift[ind] = len - 1 - j;
            }

            int ind = data(i + len - 3) * 16 + data(i + len - 2) * 4 + data(i + len - 1);
            sh1 = shift[ind];
            shift[ind] = 0;
            
            if (!sh1) sh1 = 1;

            int j = len - 1;

            while(true) {
                int sh = 1;
                while (sh && j < size) {
                    int ind = data(j - 2) * 16 + data(j - 1) * 4 + data(j);
                    sh = shift[ind];
                    j += sh;
                } 
                
                if (j < size) {
                    int is_eq = 1;
                    
                    // SIMD comparison
                    using simd_type = Kokkos::Experimental::simd<int8_t>;
                    constexpr int simd_width = simd_type::size();
                    const int cycles = len / simd_width;
                    const int leftover = len % simd_width;
                    
                    for(int k = 0; k < cycles; ++k) {
                        
                        // Load data using gather operation
                        Kokkos::Experimental::simd<int8_t> pattern1, pattern2;
                        pattern1.copy_from(&data(i + k * simd_width), Kokkos::Experimental::element_aligned_tag());
                        pattern2.copy_from(&data(j - len + 1 + k * simd_width), Kokkos::Experimental::element_aligned_tag());
                        is_eq = Kokkos::Experimental::all_of(pattern1==pattern2);
                        
                        if (!is_eq) {
                            break;
                        }
                    }
                    
                    // Handle leftover elements
                    for(int k = 0; k < leftover && is_eq; ++k) {
                        if (data(i + cycles * simd_width + k) != 
                            data(j - len + 1 + cycles * simd_width + k)) {
                            is_eq = 0;
                            break;
                        }
                    }

                    res += is_eq;
                    j += sh1;
                } else {
                    break;
                }
            }
            freq_view(i) = res;
        }
    );
Kokkos::fence();
    double stop = timer.seconds();

    // Copy results back to host
    auto freq_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), freq_view);
    for(int i = 0; i < size; ++i) {
        freq[i] = freq_host(i);
    }
    
    if (perf_collect) {
        auto coll_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), num_of_coll);
        auto comp_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), num_of_comp);
        std::cout << "number of collisions = " << coll_host(0) << "\n";
        std::cout << "number of comparasion = " << comp_host(0) << "\n";
        std::cout << freq.size() << " " << len << " " << stop << "\n";
    }
}
    Kokkos::finalize();
    return;
}