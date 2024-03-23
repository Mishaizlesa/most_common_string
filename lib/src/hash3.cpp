#include "algorithms.h"
#include "common_defs.h"
#include <immintrin.h>
typedef long long ll;
//#define COMPARASION_VECTORIZED
extern void hash3(std::vector<uint32_t>& freq ,const std::string& input_file, const uint32_t len_, const bool perf_collect) {
    std::ifstream fin(input_file);
    std::string data_;
    fin>>data_;
    ll size=data_.size();
    int len=len_;

    //vsetvlmax_e8m8 ();
    size_t VECTOR_LENGHT = 64;
    uint8_t cycles = len/VECTOR_LENGHT;
    uint8_t leftover = len%VECTOR_LENGHT;


    std::unordered_map<int8_t,int8_t>symbols_code{{'A',uint8_t(0)}, {'C',uint8_t(1)},{'G',uint8_t(2)},{'T',uint8_t(3)}};
    freq.resize(size);
    uint8_t data[size];
    for(int i=0;i<size;++i){
        data[i]=symbols_code[data_[i]];
    }
    double start = omp_get_wtime();

#pragma omp parallel for shared(freq, data) proc_bind(spread)
    for(int i=0;i<=size-len;++i){
        int res=0;
        int sh1;
        std::vector<uint32_t>shift(64,len-2);
        for(int j=2;j<len-1;++j){
            int ind=data[i+j-2]*16+data[i+j-1]*4+data[i+j];
            shift[ind]=len-1-j;
        }

        int ind=data[i+len-3]*16+data[i+len-2]*4+data[i+len-1];
        sh1=shift[ind];
        shift[ind]=0;
        
        if (!sh1) sh1=1;

        int j=len-1;

        for(;;){
            int sh=1;
            while (sh && j<size) {
                int ind=data[j-2]*16+data[j-1]*4+data[j];
                sh=shift[ind];
                j+=sh;
            } 
            if (j<size){
                int is_eq=1;
                uint8_t* pattern_1 = &data[i];
                uint8_t* pattern_2 = &data[j-len+1];
            #ifndef COMPARASION_VECTORIZED
                for(int k=0;k<len;++k){
                    if (pattern_1[k]!=pattern_2[k]){
                        is_eq=0;
                    }
                }
            #else
                for(int k=0;k<cycles;++k){
                    
                    __m512i vpattern_1 = _mm512_loadu_epi8 ((void*) pattern_1);                
                    __m512i vpattern_2 = _mm512_loadu_epi8 ((void*) pattern_2);
                    __mmask64 pcmp _mm512_cmpeq_epu8_mask (vpattern_1, vpattern_2);
                    unsigned __int64 bitmask =  _cvtmask64_u64 (pcmp);
                    if (bitmask != 0xffffffffffffffffU) {
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
                j+=sh1;

            }else{
                break;
            }
        }
        freq[i]=res;
    }
    double stop = omp_get_wtime();
    if (perf_collect) {
        std::cout<<freq.size()<<" "<<len<<" "<< stop - start << "\n";
    }
    return;
}
