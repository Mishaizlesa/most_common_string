#include "algorithms.h"
#include "common_defs.h"
//#include "riscv_vector.h"
//const size_t VECTOR_LENGHT = 64;
//#define SWAR_VECORIZED
#define COMPARASION_VECTORIZED


typedef vuint8m1_t vu8;

extern void rabin_karp_SWAR(std::vector<uint32_t>& freq ,const std::string& input_file, const uint32_t len_, const bool perf_collect) {
    std::ifstream fin(input_file);
    std::string data_;
    fin>>data_;
    uint64_t size=data_.size();
    int len=len_;

    size_t VECTOR_LENGHT = vsetvlmax_e8m1();
    size_t VECTOR_LENGHT_SWAR = vsetvlmax_e8m1();
    //std::cout<<VECTOR_LENGHT<<"\n";
    int32_t cycles = len/VECTOR_LENGHT;
    int32_t leftover = len%VECTOR_LENGHT;

   

    //std::cout << size << " " << len << std::endl;
    freq.resize(size,0);
    uint8_t data[size+VECTOR_LENGHT+len+10];
    for(int i=size;i<=size+VECTOR_LENGHT+len;++i) data[i]=5;
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
#pragma omp parallel shared(data, freq) proc_bind(close) //reduction(+: num_of_coll, num_of_comp)
{

    uint8_t eq[VECTOR_LENGHT_SWAR];

    #pragma omp for
    for (int i=0;i<size-len+1;i++){
    

        uint32_t res=0;
        uint8_t pattern_first1 = data[i];
        uint8_t pattern_first2 = data[i+1];
        uint8_t pattern_last1 = data[i+len-1];
        uint8_t pattern_last2 = data[i+len-2];
        for (int j=0;j<size-len+1;j+=VECTOR_LENGHT_SWAR)
        {

           vuint8m1_t vfirst_sym1 = vle8_v_u8m1(&data[j], VECTOR_LENGHT_SWAR);

           vuint8m1_t vfirst_sym2 = vle8_v_u8m1(&data[j+1], VECTOR_LENGHT_SWAR);

           vuint8m1_t vlast_sym1 = vle8_v_u8m1(&data[j+len-1], VECTOR_LENGHT_SWAR);

           vuint8m1_t vlast_sym2 = vle8_v_u8m1(&data[j+len-2], VECTOR_LENGHT_SWAR);

           vuint8m1_t eq_first1 =  vxor_vx_u8m1 (vfirst_sym1, pattern_first1, VECTOR_LENGHT_SWAR);

           vuint8m1_t eq_first2 =  vxor_vx_u8m1 (vfirst_sym2, pattern_first2, VECTOR_LENGHT_SWAR);

            vuint8m1_t eq_last1 =  vxor_vx_u8m1 (vlast_sym1, pattern_last1, VECTOR_LENGHT_SWAR);
            vuint8m1_t eq_last2 =  vxor_vx_u8m1 (vlast_sym2, pattern_last2, VECTOR_LENGHT_SWAR);

            vuint8m1_t veq_first = vor_vv_u8m1 (eq_first1, eq_first2, VECTOR_LENGHT_SWAR);
            vuint8m1_t veq_last = vor_vv_u8m1 (eq_last1, eq_last2, VECTOR_LENGHT_SWAR);

            vuint8m1_t veq = vor_vv_u8m1 (veq_first, veq_last, VECTOR_LENGHT_SWAR);


            vse8_v_u8m1 (eq, veq, VECTOR_LENGHT_SWAR);


            for(int t=0;t<VECTOR_LENGHT_SWAR;++t){

                if (eq[t]==0) {
                    uint8_t is_eq =1;
                    //num_of_coll++;
                #ifndef COMPARASION_VECTORIZED
                    for(int k=0;k<len;++k){
                    //num_of_comp++;
                        if (data[i+k]!=data[j+k+t]){
                            is_eq=0;
                            break;
                        }
                    }
                #else
                    uint8_t* pattern_1 = &data[i];
                    uint8_t* pattern_2 = &data[j+t];
                    for(int k=0;k<cycles;++k){
                        vuint8m1_t vpattern_1 =  vle8_v_u8m1(pattern_1, VECTOR_LENGHT);
                        vuint8m1_t vpattern_2 =  vle8_v_u8m1(pattern_2, VECTOR_LENGHT);
 
                        vbool8_t vres =  vmseq_vv_u8m1_b8(vpattern_1, vpattern_2, VECTOR_LENGHT);
 
                        uint32_t res = vmpopc_m_b8 (vres, VECTOR_LENGHT);
 
                        if (res!=VECTOR_LENGHT) {
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
#pragma omp parallel for shared(data, freq) proc_bind(close)
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
                        if (data[i+k]!=data[j+k+t]){
                            is_eq=0;
                            break;
                        }
                    }
                #else
                    uint8_t* pattern_1 = &data[i];
                    uint8_t* pattern_2 = &data[j];
                    for(int k=0;k<cycles;++k){
                        vuint8m1_t vpattern_1 =  vle8_v_u8m1(pattern_1, VECTOR_LENGHT);
                        vuint8m1_t vpattern_2 =  vle8_v_u8m1(pattern_2, VECTOR_LENGHT);
 
                        vbool8_t vres =  vmseq_vv_u8m1_b8(vpattern_1, vpattern_2, VECTOR_LENGHT);
 
                        uint32_t res = vmpopc_m_b8 (vres, VECTOR_LENGHT);
 
                        if (res!=VECTOR_LENGHT) {
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
        std::cout<<"number of collisions = "<<num_of_coll<<"\n";
        std::cout<<"number of comparasion = "<<num_of_comp<<"\n";
        std::cout<<freq.size()<<" "<<len<<" "<< stop - start << "\n";
    }
    return;
}

