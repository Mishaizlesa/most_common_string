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
extern void naive_scalar(std::vector<uint32_t>& freq, const std::string& input_file, 
                 const uint32_t len_, const bool perf_collect) {
    std::ifstream fin(input_file);
    Kokkos::Timer timer;
    int ord[256];
    std::string data;
    fin>>data;
    ll size=data.size();
    int len=len_;
    //std::cout<<size;
    len=(len<3?3:len);
    freq.resize(size);
    timer.reset();
    double st=timer.seconds();
#pragma omp parallel for shared(freq)
    for(int i=0;i<=size-len;++i){
    for(int j=0;j<size-len+1;++j){
                int is_eq=1;
                for(int k=0;k<len && is_eq;++k){
                    if (data[i+k]!=data[j+k]) is_eq=0;
                }
                freq[i]+=is_eq;
            }
        }
    std::cout<<timer.seconds()-st<<" ";
}