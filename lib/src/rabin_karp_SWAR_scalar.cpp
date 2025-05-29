#include "algorithms.h"
#include "common_defs.h"
#include <Kokkos_Core.hpp>
#include <Kokkos_SIMD.hpp>
#include <unordered_map>
#include <Kokkos_SIMD_AVX512.hpp>

extern void rabin_karp_SWAR_scalar(std::vector<uint32_t>& freq, const std::string& input_file,
                          const uint32_t len_, const bool perf_collect) {
    Kokkos::initialize();
    {
        Kokkos::Timer timer;

        std::ifstream fin(input_file);
        std::string data_;
        fin >> data_;
        const uint64_t size = data_.size();
        const uint32_t len = len_;

        // Use native SIMD width instead of trying to get it from the type
        const int vector_width = Kokkos::Experimental::basic_simd<std::uint32_t,  Kokkos::Experimental::simd_abi::avx512_fixed_size<16>>::size();
        const int swar_width = vector_width;

        // Создаем копию map для использования в device коде
        const std::unordered_map<int8_t, int8_t> host_map = {
            {'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}
        };
        
        Kokkos::View<int8_t[128]> mapSymbToCode("mapSymbToCode");
        Kokkos::parallel_for("init_map", 128, KOKKOS_LAMBDA(int i) {
            mapSymbToCode(i) = -1; // Инициализация
        });
        
        // Копируем значения из std::map в View
        for(const auto& pair : host_map) {
            mapSymbToCode(pair.first) = pair.second;
        }

        Kokkos::View<uint8_t*> data("data", size + len + vector_width);
        Kokkos::View<uint32_t*> freq_k("freq", size);

        // Инициализация данных
        Kokkos::parallel_for("init_data", size, KOKKOS_LAMBDA(int i) {
            data(i) = mapSymbToCode(data_[i]);
        });

        timer.reset();
        Kokkos::parallel_for("rabin_karp", Kokkos::RangePolicy<>(0, size-len+1),
        KOKKOS_LAMBDA(const int i) {
            uint32_t res = 0;
            uint8_t p1 = data(i);
            uint8_t p2 = data(i+1);
            uint8_t pn1 = data(i+len-1);
            uint8_t pn2 = data(i+len-2);
            
            for (uint64_t j = 0; j < size - len + 1; ++j) {
                uint64_t ptr = j;
                if (data[ptr] == p1 && 
                    data[ptr + 1] == p2 && 
                    data[ptr + len - 1] == pn1 && 
                    data[ptr + len - 2] == pn2) {
                    
                    bool is_eq = true;
                    for (uint32_t t = 0; t < len; ++t) {
                        if (data(i + t) != data(ptr + t)) {
                            is_eq = false;
                            break;
                        }
                    }
                    res += is_eq;
                }
            }
            
            freq_k(i) = res;
        });
        Kokkos::fence();
        double stop = timer.seconds();

        // Копирование результатов
        freq.resize(size);
        Kokkos::deep_copy(Kokkos::View<uint32_t*, Kokkos::HostSpace>(freq.data(), size), freq_k);

        if (perf_collect) {
            std::cout << freq.size() << " " << len << " " << stop << "\n";
        }
    }
    Kokkos::finalize();
}