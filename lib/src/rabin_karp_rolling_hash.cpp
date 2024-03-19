#include "algorithms.h"
#include "common_defs.h"

//#define HASH_VECORIZED
//#define COMPARASION_VECTORIZED
//#include "riscv_vector.h"
//const size_t VECTOR_LENGHT = 64;


extern void rabin_karp_rolling_hash(std::vector<uint32_t>& freq ,const std::string& input_file, const uint32_t len_, const bool perf_collect) {
    std::ifstream fin(input_file);
    std::string data_;
    fin>>data_;
    uint64_t size=data_.size();
    int len=len_;

    size_t VECTOR_LENGHT = 32;

    size_t VECTOR_LENGHT_HASH = 4;
    //std::cout<<VECTOR_LENGHT<<"\n";
    int32_t cycles = len/VECTOR_LENGHT;
    int32_t leftover = len%VECTOR_LENGHT;


    //std::cout << size << " " << len << std::endl;
    freq.resize(size,0);
    uint8_t data[size+len];
    for(int i=size;i<size+len;++i) data[i]=7;
    std::unordered_map<int8_t, int8_t> mapSymbToCode = { {'A', (int8_t)0}, {'C', (int8_t)1}, {'G', (int8_t)2}, {'T', (int8_t)3} };
    // prepare data
    for (int i = 0; i < size; ++i) {
        data[i] = mapSymbToCode[data_[i]];
    }
 



    double start = omp_get_wtime();

    //std::cout << timer.seconds() - st << std::endl;
#ifdef HASH_VECORIZED
    #pragma omp parallel shared(data, freq) proc_bind(close)
{

    uint64_t hash_text[VECTOR_LENGHT_HASH];


    #pragma omp for
    for (int i=0;i<size-len+1;++i){
    
        int res = 0;
        uint64_t hash_pattern = 0;
        vu64 vhash_text = vmv_v_x_u64m2 (0, VECTOR_LENGHT_HASH);
        uint64_t p = 2;
        uint64_t powmod = 1;

       for (int j = 0; j < len; ++j) {
            vu64 vdata =  vlsbu_v_u64m2 ((uint64_t*)(&data[j]), size/VECTOR_LENGHT_HASH, VECTOR_LENGHT_HASH);
            hash_pattern = hash_pattern * p + data[i + j]; 
            vhash_text =  vmadd_vx_u64m2 (vhash_text, p, vdata, VECTOR_LENGHT_HASH);
            powmod = powmod * p;
        }
        
        for (int j = 0; j < size/VECTOR_LENGHT_HASH; ++j) {
            vse64_v_u64m2 (hash_text, vhash_text, VECTOR_LENGHT_HASH);

            for(int t=0;t<VECTOR_LENGHT_HASH;++t){
                if (hash_text[t]==hash_pattern) {
                    int is_eq=1;
                    uint8_t* pattern_1 = &data[i];
                    uint8_t* pattern_2 = &data[j+(size*t)/VECTOR_LENGHT_HASH];
    #ifdef COMPARASION_VECTORIZED
                    for(int k=0;k<cycles;++k){


                        vu8 vpattern_1 =  vle8_v_u8m8(pattern_1, VECTOR_LENGHT);
                        vu8 vpattern_2 =  vle8_v_u8m8(pattern_2, VECTOR_LENGHT);


                        vbool1_t vres =  vmseq_vv_u8m8_b1(vpattern_1, vpattern_2, VECTOR_LENGHT);


                        uint32_t eq_res = vmpopc_m_b1 (vres, VECTOR_LENGHT);

            

                        if (eq_res!=VECTOR_LENGHT) {
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


            
            vu64 vdata1 =  vlsbu_v_u64m2 ((uint64_t*)(&data[j]), size/VECTOR_LENGHT_HASH, VECTOR_LENGHT_HASH);
            vu64 vdata2 =  vlsbu_v_u64m2 ((uint64_t*)(&data[j+len]), size/VECTOR_LENGHT_HASH, VECTOR_LENGHT_HASH);


            vdata1 =  vmadd_vx_u64m2 (vdata1, -powmod, vdata2, VECTOR_LENGHT_HASH);

            vhash_text =  vnmsub_vx_u64m2 (vhash_text, -p, vdata1, VECTOR_LENGHT_HASH);

        }
        freq[i] = res;
    }
}
#else
#pragma omp parallel for shared(data, freq) proc_bind(close)
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
                for(int k=0;k<cycles;++k){
                    vu8 vpattern_1 =  vle8_v_u8m8(pattern_1, VECTOR_LENGHT);
                    vu8 vpattern_2 =  vle8_v_u8m8(pattern_2, VECTOR_LENGHT);


                    vbool1_t vres =  vmseq_vv_u8m8_b1(vpattern_1, vpattern_2, VECTOR_LENGHT);


                    uint32_t eq_res = vmpopc_m_b1 (vres, VECTOR_LENGHT);

        

                    if (eq_res!=VECTOR_LENGHT) {
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

