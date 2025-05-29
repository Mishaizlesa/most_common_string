#include "algorithms.h"
#include "common_defs.h"
#include <Kokkos_Core.hpp>
#include <Kokkos_SIMD.hpp>
#include <unordered_map>

extern void rabin_karp_SWAR(std::vector<uint32_t>& freq, const std::string& input_file,
                          const uint32_t len_, const bool perf_collect) {
    Kokkos::initialize();
    {
        Kokkos::Timer timer;

    std::ifstream fin(input_file);
    std::string data_;
    fin >> data_;
    const uint64_t size = data_.size();
    const uint32_t len = len_;

    using simd_type = Kokkos::Experimental::simd<uint8_t>;
    constexpr int vector_width = simd_type::size();
    constexpr int swar_width = vector_width;

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

    int num_of_iter_search = (size - len + 1) / swar_width;
    int leftover_search = (size - len + 1) % swar_width;

    int num_of_iter_comp = (len) / swar_width;
    int leftover_comp = (len) % swar_width;

    timer.reset();
    Kokkos::parallel_for("rabin_karp", Kokkos::RangePolicy<>(0, size-len+1),
    KOKKOS_LAMBDA(const int i) {
        uint32_t res = 0;
        uint8_t p1 = data(i);
        uint8_t p2 = data(i+1);
        uint8_t pn1 = data(i+len-1);
        uint8_t pn2 = data(i+len-2);
        int f = 1;
        for(uint64_t j = 0, it = 0; it < num_of_iter_search; j += swar_width, it++) {
            
            simd_type v1, v2, vn1, vn2;
                
            v1.copy_from(&data(j), Kokkos::Experimental::element_aligned_tag());
            v2.copy_from(&data(j + 1), Kokkos::Experimental::element_aligned_tag());
            vn1.copy_from(&data(j + len - 1), Kokkos::Experimental::element_aligned_tag());
            vn2.copy_from(&data(j + len - 2), Kokkos::Experimental::element_aligned_tag());
            auto mask = (v1 == p1) && (v2 == p2) && (vn1 == pn1) && (vn2 == pn2); 
            
            // Проверка полного совпадения
                for(uint64_t k = 0; k < swar_width; ++k) {
                    bool is_eq = 0;
                    if (mask[k]) {
                        is_eq = 1;
                        for (uint32_t it2 = 0, t = 0; it2 < num_of_iter_comp && is_eq; it2++, t += swar_width){
                            v1.copy_from(&data(i + swar_width * t), Kokkos::Experimental::element_aligned_tag());
                            v2.copy_from(&data(j + k + swar_width * t), Kokkos::Experimental::element_aligned_tag());
                            is_eq = Kokkos::Experimental::all_of(v1==v2);
                        }
                        for (int t = 0; t < leftover_comp &&  is_eq; ++t)
                        {
                            if (data(i + swar_width * num_of_iter_comp + t) != data(j + k + swar_width * num_of_iter_comp + t)){
                                is_eq = 0;
                                break;
                            }
                        }
                    }
                    res += is_eq;
                }
        }
        for (int j = 0; j < leftover_search; ++j)
        {
            int ptr = num_of_iter_search * swar_width + j;
            if (data[ptr] == p1 && data[ptr + 1] == p2 && data[ptr + len - 1] == pn1 && data[ptr + len - 2] == pn2)
            {
                simd_type v1, v2;
                bool is_eq = 1;
                is_eq = 1;
                for (uint32_t it2 = 0, t = 0; it2 < num_of_iter_comp && is_eq; it2++, t += swar_width){
                    v1.copy_from(&data(i + swar_width * t), Kokkos::Experimental::element_aligned_tag());
                    v2.copy_from(&data(ptr + swar_width * t), Kokkos::Experimental::element_aligned_tag());
                    is_eq = Kokkos::Experimental::all_of(v1==v2);
                }
                for (int t = 0; t < leftover_comp; ++t)
                {
                    if (data(i + swar_width * num_of_iter_comp + t) != data(ptr + swar_width * num_of_iter_comp + t)){
                        is_eq = 0;
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