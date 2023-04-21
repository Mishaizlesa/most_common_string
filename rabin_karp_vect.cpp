#include <Kokkos_Core.hpp>
#include <Kokkos_Atomic.hpp>
#include <omp.h>
#include <cstdio>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
typedef long long ll;
ll mod=1e9+7;
ll p=31;
int main(int argc, char* argv[]){
    std::ifstream fin(argv[1]);
    int f=argv[2][0]-'0';
    std::ofstream fout("tmp.txt");
    Kokkos::initialize(argc, argv);{
        Kokkos::View<ll*, Kokkos::SharedSpace> ord("ord", 256);
        ord['A']=1LL;
        ord['C']=2LL;
        ord['G']=3LL;
        ord['T']=4LL;
        std::string data_;
        fin>>data_;
        ll size=data_.size();
        int len=(data_.size()>10000?1000:data_.size()/10);
        len=(len<3?3:len);
        Kokkos::View<double*,Kokkos::SharedSpace> freq("array", size);
        Kokkos::View<char*,Kokkos::SharedSpace> data("device_string",size);
        for(int i=0;i<size;++i){
            data[i]=data_[i];
        }
        ll def_val=0;
        ll powmod=1;
        for(int i=0;i<len;++i){
            def_val=(def_val+ord[data[i]]*powmod)%mod;
            powmod*=p;
            powmod%=mod;
        }
        Kokkos::Timer timer;
        timer.reset();
        double st = timer.seconds();
        Kokkos::parallel_for( "yAx", size-len+1, KOKKOS_LAMBDA (int i) {
            ll res=0;
            ll hash=0;
            ll p=31;
            ll mod=1e9+7;
            ll powmod=1;
            ll tmp_h=0;
            for(int j=0;j<len;++j){
                hash=(hash+ord[data[i+j]]*powmod)%mod;
                powmod*=p;
                powmod%=mod;
            }
            tmp_h=def_val;
            if (tmp_h==hash){
                int is_eq=1;
                for(int k=0;k<len;++k){
                    if (data[i+k]!=data[k]){
                        is_eq=0;break;
                    }
                }
                res+=is_eq;
            }
            ll powlen=1;
            powmod=1;
            for(int k=0;k<len;++k){
                powlen*=p;
                powlen%=mod;
            }
            for(int j=1;j<=size-len;++j){
                tmp_h=(tmp_h-(ord[data[j-1]]*powmod)%mod+(ord[data[len+j-1]]*powmod*powlen)%mod+mod)%mod;
                powmod*=p;
                powmod%=mod;
                if (tmp_h==(hash*powmod)%mod){
                    int is_eq=1;
                    for(int k=0;k<len;++k){
                        if (data[i+k]!=data[j+k]){
                            is_eq=0;break;
                        }
                    }
                    res+=is_eq;
                }
            }
            freq(i)=res;
        });
        if (f==1) fout<<timer.seconds()-st<<" ";
        else if (f==2){
            for(int i=0;i<size-len+1;++i){
                fout<<freq[i];
            }
        }
    }
    Kokkos::finalize();
    return 0;
}
