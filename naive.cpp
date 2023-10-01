#include <omp.h>
#include <iostream>
#include <string>
#include <fstream>
#include <Kokkos_Core.hpp>
typedef long long ll;
int main(int argc, char* argv[]){
    std::ifstream fin(argv[1]);
    int f=argv[2][0]-'0';
    std::ofstream fout("tmp.txt");
    Kokkos::Timer timer;
    timer.reset();
    Kokkos::initialize(argc, argv);{
        std::string data_;
        fin>>data_;
        ll size=data_.size();
        int len=std::atoi(argv[3]);
        Kokkos::View<double*, Kokkos::SharedSpace> freq("array", size);
        Kokkos::View<char*,Kokkos::SharedSpace> data("device_string",size);
        for(int i=0;i<size;++i){
            data[i]=data_[i];
        }
        Kokkos::parallel_for( "yAx", size-len+1, KOKKOS_LAMBDA (int i) {
            for(int j=0;j<size-len+1;++j){
                int is_eq=1;
                for(int k=0;k<len && is_eq;++k){
                    if (data[i+k]!=data[j+k]) is_eq=0;
                }
                freq[i]+=is_eq;
            }
        });
	Kokkos::fence();
        if (f==1){
            fout<<timer.seconds()<<" ";
        }
        if (f==2){
            std::ofstream fout(argv[3]);
            for(int i=0;i<size-len+1;++i){
                fout<<freq[i];
            }
        }
    }
    Kokkos::finalize();
    return 0;
}
