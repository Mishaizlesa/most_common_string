#include "algorithms.h"
#include "common_defs.h"
//#include "riscv_vector.h"
//const size_t VECTOR_LENGHT = 64;
extern void rabin_karp(std::vector<uint32_t>& freq ,const std::string& input_file, const uint32_t len_, const bool perf_collect) {
    std::ifstream fin(input_file);
    std::string data_;
    fin>>data_;
    uint64_t size=data_.size();
    int len=len_;

    size_t VECTOR_LENGHT = vsetvlmax_e8m1();
    //std::cout<<VECTOR_LENGHT<<"\n";
    int32_t cycles = len/VECTOR_LENGHT;
    int32_t rest = len%VECTOR_LENGHT;


    //std::cout << size << " " << len << std::endl;
    freq.resize(size,0);
    int8_t data[size+1];
    data[size] = 'A';
    std::unordered_map<int8_t, int8_t> mapSymbToCode = { {'A', (int8_t)0}, {'C', (int8_t)1}, {'G', (int8_t)2}, {'T', (int8_t)3} };
    // prepare data
    for (int i = 0; i < size; ++i) {
        data[i] = mapSymbToCode[data_[i]];
    }



    double start = omp_get_wtime();
    //std::cout << timer.seconds() - st << std::endl;
#pragma omp parallel for shared(data, freq)
    for (int i=0;i<size-len+1;++i){
        //for (int i = 0; i < size - len + 1; i++) {
        int res = 0;
        uint64_t hash_pattern = 0;
        uint64_t hash_text = 0;
        uint64_t mod = (1LL<<31) - 1;
        uint64_t p = 2;
        uint64_t powmod = 1;
        for (int j = 0; j < len; ++j) {
            hash_pattern = hash_pattern * p + data[i + j];
            hash_text = hash_text * p + data[j];
            powmod = powmod * p;
            //powmod%=mod;
            //hash_pattern %= mod;
            //hash_text %= mod;
        }
        for (int j = 0; j < size - len + 1; ++j) {
            if (hash_text == hash_pattern) {
                int is_eq=1;
                int8_t* pattern_1 = &data[i];
                int8_t* pattern_2 = &data[j];
                for(int k=0;k<cycles;++k){
                    vint8m1_t vpattern_1 =  vle8_v_i8m1(pattern_1, VECTOR_LENGHT);
                    vint8m1_t vpattern_2 =  vle8_v_i8m1(pattern_2, VECTOR_LENGHT);

                    vbool8_t vres =  vmseq_vv_i8m1_b8(vpattern_1, vpattern_2, VECTOR_LENGHT);

                    uint32_t res = vmpopc_m_b8 (vres, VECTOR_LENGHT);

                    if (res!=VECTOR_LENGHT) {
                        is_eq=0;
                        break;
                    }

                    pattern_1+=VECTOR_LENGHT;
                    pattern_2+=VECTOR_LENGHT;


                }
                /*VECTOR_LENGHT = vsetvl_e8m8(rest);
                //printf("%d %d\n", VECTOR_LENGHT, rest);
                vint8m8_t vpattern_1 =  vle8_v_i8m8(pattern_1, VECTOR_LENGHT);
                vint8m8_t vpattern_2 =  vle8_v_i8m8(pattern_2, VECTOR_LENGHT);

                vbool1_t vres =  vmseq_vv_i8m8_b1(vpattern_1, vpattern_2, VECTOR_LENGHT);

                int32_t num_of_eq = vmpopc_m_b1 (vres, VECTOR_LENGHT);
                //printf("%d %d %d %d\n", res, rest, is_eq, res==rest);
                if (num_of_eq!=VECTOR_LENGHT) {
                    is_eq=0;
                }*/
                for(int k=0;k<rest && is_eq;++k){
                    if (*pattern_1!=*pattern_2){
                        is_eq=0;
                        break;
                    }
                    pattern_1++;
                    pattern_2++;
                }
                /*for(int k=0;k<len && is_eq;++k){
                    if (data[i+k]!=data[j+k]) is_eq=0;
                }*/
                res+=is_eq;
            }
            hash_text = (hash_text * p - data[j] * powmod + data[j + len]);
        }
        freq[i] = res;
    }

    double stop = omp_get_wtime();
    if (perf_collect) {
        std::cout<<freq.size()<<" "<<len<<" "<< stop - start << "\n";
    }
    return;
}