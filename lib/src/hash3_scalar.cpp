#include <Kokkos_Core.hpp>
#include <Kokkos_Atomic.hpp>
#include <omp.h>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <fstream>
typedef long long ll;
extern void hash3_scalar(std::vector<uint32_t>& freq, const std::string& input_file, 
                 const uint32_t len_, const bool perf_collect) {
    std::ifstream fin(input_file);
    Kokkos::Timer timer;
    int ord[256];
    std::string data;
    fin>>data;
    ll size=data.size();
    int len=len_;
    std::vector<uint8_t> data_(size);
    //std::cout<<size;
    len=(len<3?3:len);
    ord['A']=0;
    ord['C']=1;
    ord['G']=2;
    ord['T']=3;
    for (int i=0; i < size; ++i)
    {
        data_[i] = ord[data[i]];
    }
    freq.resize(size);
    timer.reset();
    double st=timer.seconds();
#pragma omp parallel for shared(freq)
    for(int i=0;i<=size-len;++i){
        int res=0;
        int sh1;
        std::vector<int>shift(64,len-2);
        ll hash=0;
        for(int j=2;j<=len-1;++j){
            int ind=data_[i+j-2]*16+data_[i+j-1]*4+data_[i+j];
            if (j==len-1) sh1=shift[ind];
            shift[ind]=len-1-j;
        }
        
        if (!sh1) sh1=1;
        int j=len-1;
        for(;;){
            int sh=1;
            while (sh && j<size) {
                int ind=data_[j-2]*16+data_[j-1]*4+data_[j];
                sh=shift[ind];
                j+=sh;
            }
            if (j<size){
                int is_eq=1;
                for(int k=0;k<len;++k){
                    if (data_[i+k]!=data_[j-len+1+k]){
                        is_eq=0;
                        break;
                    }
                }
                res+=is_eq;
                j+=sh1;
            }else{
                break;
            }
        }
        freq[i]=res;
    }
    std::cout<<timer.seconds()-st<<" ";
}