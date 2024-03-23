#include "algorithms.h"
#include "common_defs.h"
#include <immintrin.h>

//#define HASH_VECORIZED
#define COMPARASION_VECTORIZED

extern void rabin_karp_rolling_hash(std::vector<uint32_t>& freq ,const std::string& input_file, const uint32_t len_, const bool perf_collect) {
    std::ifstream fin(input_file);
    std::string data_;
    fin>>data_;
    uint64_t size=data_.size();
    int len=len_;

    size_t VECTOR_LENGHT = 64;

    size_t VECTOR_LENGHT_HASH = 8;
    //std::cout<<VECTOR_LENGHT<<"\n";
    int32_t cycles = len/VECTOR_LENGHT;
    int32_t leftover = len%VECTOR_LENGHT;


    freq.resize(size,0);
    uint8_t data[size+8-size%8];
    uint64_t data_load[size+8-size%8];
    for(int i=size;i<size+8-size%8;++i){
	data[i]=7;
	data_load[i]=7;
   }
    std::unordered_map<int8_t, int8_t> mapSymbToCode = { {'A', (int8_t)0}, {'C', (int8_t)1}, {'G', (int8_t)2}, {'T', (int8_t)3} };
    // prepare data
    for (int i = 0; i < size; ++i) {
        data[i] = mapSymbToCode[data_[i]];
	data_load[i]=data[i];
    }
 

    uint32_t indicies[8];
    uint64_t size_ = size+8-size%8;
    for(int i=0;i<VECTOR_LENGHT_HASH;++i){
		
	indicies[i] = (i*size_);

    }
   
    double start = omp_get_wtime();

    //std::cout << timer.seconds() - st << std::endl;
#ifdef HASH_VECORIZED
    #pragma omp parallel shared(data, freq) proc_bind(spread)
{

    uint64_t hash_text[VECTOR_LENGHT_HASH];
    uint64_t p = 2;
    __m256i vindices = _mm256_loadu_epi32 ((void*) indicies);   
    __m512i vp = _mm512_set1_epi64(p);

    #pragma omp for
    for (int i=0;i<size-len+1;++i){
    
        int res = 0;
        uint64_t hash_pattern = 0;
        __m512i vhash_text = _mm512_set1_epi64(0);
        uint64_t powmod = 1;

       for (int j = 0; j < len; ++j) {
            __m512i vdata =  _mm512_i32gather_epi64(vindices, &data_load[j], 1);
            hash_pattern = hash_pattern * p + data[i + j]; 
            vhash_text =  _mm512_mullo_epi64(vhash_text, vp);
	    vhash_text = _mm512_add_epi64 (vhash_text, vdata);
            powmod = powmod * p;
        }
	__m512i vhash_pattern = _mm512_set1_epi64(hash_pattern);
	__m512i vpowmod = _mm512_set1_epi64(powmod);

        
        for (int j = 0; j < (size_)/VECTOR_LENGHT_HASH; ++j) {
             __mmask8 pcmp =  _mm512_cmpeq_epu64_mask  (vhash_text, vhash_pattern);
            unsigned int bitmask =  _cvtmask8_u32 (pcmp);  
            for(int t=0;t<VECTOR_LENGHT_HASH;++t){
                if (bitmask & (1<<(t))) {
                    int is_eq=1;
                    uint8_t* pattern_1 = &data[i];
                    uint8_t* pattern_2 = &data[j+(size_*t)/VECTOR_LENGHT_HASH];
    #ifdef COMPARASION_VECTORIZED
                 /*for(int k=0;k<cycles;++k){
                    
                    __m512i vpattern_1 =  _mm512_loadu_epi8 ((void*) pattern_1);            
                    __m512i vpattern_2 =  _mm512_loadu_epi8 ((void*) pattern_2);
                    __mmask64 pcmp =  _mm512_cmpeq_epu8_mask  (vpattern_1, vpattern_2);
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
                }*/
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

    #else
                    for(int k=0;k<len;++k){
                        if (pattern_1[k]!=pattern_2[k]){
                            is_eq=0;
                            break;
                        }
                    }
    #endif
                    res+=is_eq;
                }
            }


            
            __m512i vdata1 =  _mm512_i32gather_epi64(vindices, (void*) &data_load[j], 1);
            __m512i vdata2 =  _mm512_i32gather_epi64(vindices, (void*) &data_load[j+len], 1);



            vdata1 =  _mm512_mullo_epi64(vdata1, vpowmod);
	    vdata1 = _mm512_sub_epi64 (vdata1, vdata2);


            vhash_text =  _mm512_mullo_epi64(vhash_text, vp);
	    vhash_text =  _mm512_sub_epi64 (vhash_text, vdata1);

        }
        freq[i] = res;
    }
}
#else
#pragma omp parallel for shared(data, freq) proc_bind(spread)
    for (int i=0;i<size-len+1;++i){
        int res = 0;
        uint64_t hash_pattern = 0;
        uint64_t hash_text = 0;
        uint64_t p = 2;
        uint64_t powmod = 1;
        for (int j = 0; j < len; ++j) {
            hash_pattern = hash_pattern * p + data[i + j];
            hash_text = hash_text * p + data[j];
            powmod = powmod * p;
        }
        for (int j = 0; j < size - len + 1; ++j) {
            if (hash_text == hash_pattern) {
                int is_eq=1;
                uint8_t* pattern_1 = &data[i];
                uint8_t* pattern_2 = &data[j];
            #ifdef COMPARASION_VECTORIZED
                /*for(int k=0;k<cycles;++k){
                    
                    __m512i vpattern_1 =  _mm512_loadu_epi8 ((void*) pattern_1);            
                    __m512i vpattern_2 =  _mm512_loadu_epi8 ((void*) pattern_2);
                    __mmask64 pcmp =  _mm512_cmpeq_epu8_mask  (vpattern_1, vpattern_2);
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
                }*/
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


            #else
                for(int k=0;k<len;++k){
                    if (pattern_1[k]!=pattern_2[k]){
                        is_eq=0;
                        break;
                    }
                }
            #endif
                res+=is_eq;
            }
            hash_text = (hash_text * p - data[j] * powmod + data[j + len]);
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

