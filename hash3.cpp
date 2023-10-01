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
    int f=argv[2][0]-'0';
    std::ofstream fout("tmp.txt");
    Kokkos::initialize(argc, argv);{
        Kokkos::Timer timer;
        //int ord_[256];
        std::string data_;
        fin>>data_;
        ll size=data_.size();
        int len=std::atoi(argv[3]);
        len=(len<3?3:len);
        Kokkos::View<int*,Kokkos::SharedSpace> freq("array", size);
        Kokkos::View<char*,Kokkos::SharedSpace> tmp("tmp",size);
        for(int i=0;i<size;++i) tmp[i]=data_[i];
        Kokkos::View<const char*,Kokkos::SharedSpace> data=tmp;
        timer.reset();
        Kokkos::parallel_for( "yAx", size-len+1, KOKKOS_LAMBDA (int i) {
            int shift[400];
            for(int j=0;j<400;++j) shift[j]=len-2;
            int res=0;
            int sh1;
            ll hash=0;
            for(int j=2;j<=len-1;++j){
                int ind=(data[i+j-2]-'A')*16+(data[i+j-1]-'A')*4+(data[i+j]-'A');
                if (j==len-1) sh1=shift[ind];
                shift[ind]=len-1-j;
            }
            
            if (!sh1) sh1=1;
            int j=len-1;
            for(;;){
                int sh=1;
                while (sh && j<size) {
                    int ind=(data[j-2]-'A')*16+(data[j-1]-'A')*4+(data[j]-'A');
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
            freq(i)=res;
        });
        Kokkos::fence();
        if (f==1){ fout<<timer.seconds()<<" ";}
        if (f==2){
            for(int i=0;i<=size-len;++i){
                fout<<freq[i];
            }
        }else if (f==3){
            int mmax=0;
            for(int i=0;i<size;++i){
                mmax=(mmax>freq[i]?mmax:freq[i]);
            }
            fout<<mmax;
        }
    }
    Kokkos::finalize();
    return 0;
}

