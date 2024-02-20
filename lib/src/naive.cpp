#include "algorithms.h"
#include "common_defs.h"
//#include "riscv_vector.h"
typedef long long ll;
extern void naive(std::vector<uint32_t>& freq ,const std::string& input_file, const uint32_t len_, const bool perf_collect){
    std::ifstream fin(input_file);
    std::string data_;
    fin>>data_;
    ll size=data_.size();
    int len=len_;
    size_t VECTOR_LENGHT = vsetvlmax_e8m1();
    uint8_t cycles = len/VECTOR_LENGHT;
    uint8_t rest = len%VECTOR_LENGHT;


    std::unordered_map<int8_t,int8_t>symbols_code{{'A',int8_t(0)}, {'C',int8_t(1)},{'G',int8_t(2)},{'T',int8_t(3)}};
    freq.resize(size);
    int8_t data[size];
    for(int i=0;i<size;++i){
        data[i]=symbols_code[data_[i]];
    }


    double start = omp_get_wtime();
    #pragma omp parallel for shared(data, freq) proc_bind(spread)
        for (int i=0;i<size-len+1;++i){
            uint32_t res=0;
            for(int j=0;j<size-len+1;++j){
                int is_eq=1;
                //vbool1_t vmseq_vv_i8m8_b1 (vint8m8_t op1, vint8m8_t op2, size_tvl);
                //vint8m1_t vle8_v_i8m1 (const int8_t *base, size_t vl);
                //vuint8m1_t vredand_vs_u8m8_u8m1_m (vbool1_t mask, vuint8m1_tdest, vuint8m8_t vector, vuint8m1_t scalar, size_t vl);
                //unsigned long vpopc_m_b1 (vbool1_t op1, size_t vl);
                //size_t vsetvlmax_e8m1 ();
                //size_t vsetvlmax_e8m2 ();
                //size_t vsetvlmax_e8m4 ();
                //size_t vsetvlmax_e8m8 ();
                int8_t* pattern_1 = &data[i];
                int8_t* pattern_2 = &data[j];
                for(int k=0;k<cycles;++k){
                    vint8m1_t vpattern_1 =  vle8_v_i8m1(pattern_1, VECTOR_LENGHT);
                    vint8m1_t vpattern_2 =  vle8_v_i8m1(pattern_2, VECTOR_LENGHT);

                    vbool8_t vres =  vmseq_vv_i8m1_b8(vpattern_1, vpattern_2, VECTOR_LENGHT);

                    uint32_t num_of_eq = vmpopc_m_b8 (vres, VECTOR_LENGHT);

                    if (num_of_eq!=VECTOR_LENGHT) {
                        is_eq=0;
                        break;
                    }

                    pattern_1+=VECTOR_LENGHT;
                    pattern_2+=VECTOR_LENGHT;

                }
                for(int k=0;k<rest && is_eq;++k){
                    if (*pattern_1!=*pattern_2){
                        is_eq=0;
                        break;
                    }
                    pattern_1++;
                    pattern_2++;
                }
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


