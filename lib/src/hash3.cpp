#include "algorithms.h"
#include "common_defs.h"
typedef long long ll;
extern void hash3(std::vector<uint32_t>& freq ,const std::string& input_file, const uint32_t len_, const bool perf_collect) {
    std::ifstream fin(input_file);
    std::string data_;
    fin>>data_;
    ll size=data_.size();
    int len=len_;

    //vsetvlmax_e8m8 ();
    size_t VECTOR_LENGHT = vsetvlmax_e8m8();
    uint8_t cycles = len/VECTOR_LENGHT;
    uint8_t rest = len%VECTOR_LENGHT;


    std::unordered_map<int8_t,int8_t>symbols_code{{'A',uint8_t(0)}, {'C',uint8_t(1)},{'G',uint8_t(2)},{'T',uint8_t(3)}};
    freq.resize(size);
    int8_t data[size];
    for(int i=0;i<size;++i){
        data[i]=symbols_code[data_[i]];
    }
    double start = omp_get_wtime();
#pragma omp parallel for shared(freq, data) proc_bind(spread)
    for(int i=0;i<=size-len;++i){
        int res=0;
        int sh1;
        std::vector<int>shift(64,len-2);
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
                //vbool1_t vmseq_vv_i8m8_b1 (vint8m8_t op1, vint8m8_t op2, size_tvl);
                //vint8m1_t vle8_v_i8m1 (const int8_t *base, size_t vl);
                //vuint8m1_t vredand_vs_u8m8_u8m1_m (vbool1_t mask, vuint8m1_tdest, vuint8m8_t vector, vuint8m1_t scalar, size_t vl);
                //unsigned long vpopc_m_b1 (vbool1_t op1, size_t vl);
                //size_t vsetvlmax_e8m1 ();
                //size_t vsetvlmax_e8m2 ();
                //size_t vsetvlmax_e8m4 ();
                //size_t vsetvlmax_e8m8 ();
                int8_t* pattern_1 = &data[i];
                int8_t* pattern_2 = &data[j-len+1];
                for(int k=0;k<cycles;++k){
                    vint8m8_t vpattern_1 =  vle8_v_i8m8(pattern_1, VECTOR_LENGHT);
                    vint8m8_t vpattern_2 =  vle8_v_i8m8(pattern_2, VECTOR_LENGHT);

                    vbool1_t vres =  vmseq_vv_i8m8_b1(vpattern_1, vpattern_2, VECTOR_LENGHT);

                    uint32_t num_of_eq = vmpopc_m_b1 (vres, VECTOR_LENGHT);

                    if (num_of_eq!=VECTOR_LENGHT) {
                        is_eq=0;
                        break;
                    }

                    pattern_1+=VECTOR_LENGHT;
                    pattern_2+=VECTOR_LENGHT;


                }
                /*for(int k=0;k<rest && is_eq;++k){
                    if (*pattern_1!=*pattern_2){
                        is_eq=0;
                        break;
                    }
                    pattern_1++;
                    pattern_2++;
                }
                res+=is_eq;*/
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
