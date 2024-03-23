#include "algorithms.h"
#include "common_defs.h"
//#include "riscv_vector.h"
//const size_t VECTOR_LENGHT = 64;
#define SWAR_VECORIZED
#define COMPARASION_VECTORIZED



extern void rabin_karp_SWAR(std::vector<uint32_t>& freq ,const std::string& input_file, const uint32_t len_, const bool perf_collect) {
    std::ifstream fin(input_file);
    std::string data_;
    fin>>data_;
    uint64_t size=data_.size();
    int len=len_;

    size_t VECTOR_LENGHT = 32;
    size_t VECTOR_LENGHT_SWAR = 64;
    //std::cout<<VECTOR_LENGHT<<"\n";
    int32_t cycles = len/VECTOR_LENGHT;
    int32_t leftover = len%VECTOR_LENGHT;

   

    //std::cout << size << " " << len << std::endl;
    freq.resize(size,0);
    uint8_t data[size+1+64];
    for(int i=0;i<=64;++i) data[size+i] = 5;
    std::unordered_map<int8_t, int8_t> mapSymbToCode = { {'A', (int8_t)0}, {'C', (int8_t)1}, {'G', (int8_t)2}, {'T', (int8_t)3} };
    // prepare data
    for (int i = 0; i < size; ++i) {
        data[i] = mapSymbToCode[data_[i]];
    }
 



    double start = omp_get_wtime();

    uint32_t num_of_coll=0;
    uint32_t num_of_comp=0;
    //std::cout << timer.seconds() - st << std::endl;
#ifdef SWAR_VECORIZED
#pragma omp parallel shared(data, freq) proc_bind(spread)
{

    __m512i vzeroes = _mm512_set1_epi8(0);

    #pragma omp for
    for (int i=0;i<size-len+1;i++){
    

        uint32_t res=0;
	__m512i vpattern_first1 = _mm512_set1_epi8(data[i]);
	__m512i vpattern_first2 = _mm512_set1_epi8(data[i+1]);

	__m512i vpattern_last1 = _mm512_set1_epi8(data[i+len-1]);

	__m512i vpattern_last2 = _mm512_set1_epi8(data[i+len-2]);


        for (int j=0;j<size-len+1;j+=VECTOR_LENGHT_SWAR)
        {

	   __m512i vfirst_sym1 =  _mm512_loadu_epi8 ((void*) &data[j]);
	   __m512i vfirst_sym2 =  _mm512_loadu_epi8 ((void*) &data[j+1]);

	   __m512i vlast_sym1 =  _mm512_loadu_epi8 ((void*) &data[j+len-1]);
           __m512i vlast_sym2 =  _mm512_loadu_epi8 ((void*) &data[j+len-2]);
	   
 
           __m512i eq_first1 =  _mm512_xor_si512 (vfirst_sym1, vpattern_first1);
           __m512i eq_first2 =  _mm512_xor_si512 (vfirst_sym2, vpattern_first2);
           __m512i eq_last1 =  _mm512_xor_si512 (vlast_sym1, vpattern_last1);
           __m512i eq_last2 =  _mm512_xor_si512 (vlast_sym2, vpattern_last2);
           __m512i veq_first = _mm512_or_si512 (eq_first1, eq_first2);
           __m512i veq_last =  _mm512_or_si512 (eq_last1, eq_last2);
           __m512i veq = _mm512_or_si512 (veq_first, veq_last);


            __mmask64 pcmp =  _mm512_cmpeq_epi8_mask  (veq, vzeroes);
            uint64_t bitmask =  _cvtmask64_u64 (pcmp);  


            for(int t=0;t<VECTOR_LENGHT_SWAR;++t){

                if (bitmask & (1LL<<t)) {
		    uint8_t* pattern_1 = &data[i];
		    uint8_t* pattern_2 = &data[j+t];
                    uint8_t is_eq =1;
                #ifndef COMPARASION_VECTORIZED
                    for(int k=0;k<len;++k){
                        if (pattern_1[k]!=pattern_2[k]){
                            is_eq=0;
                            break;
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
            }
        }
        freq[i] = res;
    }
}
#else
#pragma omp parallel for shared(data, freq) proc_bind(spread)
    for (int i=0;i<size-len+1;i++){
        uint32_t res=0;
        uint8_t pattern_first1 = data[i];
        uint8_t pattern_first2 = data[i+1];
        uint8_t pattern_last1 = data[i+len-1];
        uint8_t pattern_last2 = data[i+len-2];
        for (int j=0;j<size-len+1;j++) {
            if (data[j]==pattern_first1 && data[j+1]==pattern_first2 && data[j+len-1]==pattern_last1 && data[j+len-2]==pattern_last2){
                uint8_t is_eq=1;
                #ifndef COMPARASION_VECTORIZED
                    for(int k=0;k<len;++k){
                    //num_of_comp++;
                        if (data[i+k]!=data[j+k]){
                            is_eq=0;
                            break;
                        }
                    }
                #else
                    uint8_t* pattern_1 = &data[i];
                    uint8_t* pattern_2 = &data[j];
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
        }
        freq[i] = res;
    }
#endif
    double stop = omp_get_wtime();
    if (perf_collect) {
        std::cout<<freq.size()<<" "<<len<<" "<< stop - start << "\n";
    }
    return;
}

