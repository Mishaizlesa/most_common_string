#include <Kokkos_Core.hpp>
#include <Kokkos_Atomic.hpp>
#include <omp.h>
#include <cstdio>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
typedef long long ll;
int main(int argc, char* argv[]) {
    std::ifstream fin(argv[1]);
    int f=(argc>2);
    std::ofstream fout("tmp.txt");
    Kokkos::initialize(argc, argv);{
        Kokkos::Timer timer;
        Kokkos::View<ll*, Kokkos::SharedSpace> ord("ord", 256);
        std::string data_;
        fin>>data_;
        ll size=data_.size();
        int len=(data_.size()>10000?1000:data_.size()/10);
        len=(len<3?3:len);
        ord['A']=0LL;
        ord['C']=1LL;
        ord['G']=2LL;
        ord['T']=3LL;
        Kokkos::View<double*, Kokkos::SharedSpace> freq("array", size);
        Kokkos::View<char*,Kokkos::SharedSpace> data("device_string",size);
        for(int i=0;i<size;++i){
            data[i]=data_[i];
        }
        timer.reset();
        double st=timer.seconds();
        Kokkos::parallel_for( "yAx", size-len+1, KOKKOS_LAMBDA (int i) {
            //std::cout<<i<<"\n";
            int res=0;
            int sh1;
            std::vector<int>shift(64,len-2);
            ll hash=0;
            for(int j=2;j<=len-1;++j){
                int ind=ord[data[i+j-2]]*16+ord[data[i+j-1]]*4+ord[data[i+j]];
                if (j==len-1) sh1=shift[ind];
                shift[ind]=len-1-j;
            }
            
            if (!sh1) sh1=1;
            int j=len-1;
            for(;;){
                int sh=1;
                while (sh && j<size) {
                    int ind=ord[data[j-2]]*16+ord[data[j-1]]*4+ord[data[j]];
                    sh=shift[ind];
                    j+=sh;
                }
                if (j<size){
                    int is_eq=1;
                    for(int k=0;k<len;++k){
                        if (data[i+k]!=data[j-len+1+k]){
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
        });
        if (f){ fout<<timer.seconds()-st<<" ";}
    }
    Kokkos::finalize();
    return 122;
}

