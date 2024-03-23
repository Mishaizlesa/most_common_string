#include "algorithms.h"
#include "common_defs.h"
#include <immintrin.h>
//#include "riscv_vector.h"
#define COMPARASION_VECTORIZED
typedef long long ll;
extern void naive(std::vector<uint32_t>& freq ,const std::string& input_file, const uint32_t len_, const bool perf_collect){
    std::ifstream fin(input_file);
    std::string data_;
    fin>>data_;
    ll size=data_.size();
    int len=len_;
    size_t VECTOR_LENGHT = 16;
    uint8_t cycles = len/VECTOR_LENGHT;
    uint8_t leftover = len%VECTOR_LENGHT;


    std::unordered_map<int8_t,int8_t>symbols_code{{'A',int8_t(0)}, {'C',int8_t(1)},{'G',int8_t(2)},{'T',int8_t(3)}};
    freq.resize(size);
    uint8_t data[size];
    for(int i=0;i<size;++i){
        data[i]=symbols_code[data_[i]];
    }


    double start = omp_get_wtime();
    #pragma omp parallel for shared(data, freq) proc_bind(spread)
        for (int i=0;i<size-len+1;++i){
            uint32_t res=0;
            for(int j=0;j<size-len+1;++j){
                int is_eq=1;
                uint8_t* pattern_1 = &data[i];
                uint8_t* pattern_2 = &data[j];
            #ifndef COMPARASION_VECTORIZED
                for(int k=0;k<len;++k){
                    if (pattern_1[k]!=pattern_2[k]){
                        is_eq=0;
                    }
                }
            #else
                for(int k=0;k<cycles;++k){
                    
                    __m128i vpattern_1 =  _mm_loadu_epi8 ((void*) pattern_1);                
                    __m128i vpattern_2 =  _mm_loadu_epi8 ((void*) pattern_2);
                    __mmask16 pcmp =  _mm_cmpeq_epu8_mask  (vpattern_1, vpattern_2);
                    unsigned bitmask =  _cvtmask16_u32 (pcmp);
                    if (bitmask != 0xffffU) {
                        is_eq=0;
                        break;
                    }

                    pattern_1+=VECTOR_LENGHT;
                    pattern_2+=VECTOR_LENGHT;

                }
                for(int k=0;k<leftover && is_eq;++k){
                    if (*pattern_1!=*pattern_2){
                        is_eq=0;
                        break;
                    }
                    pattern_1++;
                    pattern_2++;
                }
            #endif
                res+=is_eq;
            }
            freq[i]=res;
        }
    
    

    double stop = omp_get_wtime();
    if (perf_collect) {
        std::cout<<freq.size()<<" "<<len<<" "<< stop - start << "\n";
    }
    return;
}

