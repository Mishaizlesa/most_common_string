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
std::vector<ll>pre(int len){
    std::vector<ll>powmod(len);
    powmod[0]=1;
    for(int i=1;i<len;++i){
        powmod[i]=(powmod[i-1]*p)%mod;
    }
    return powmod;
}
int main(int argc, char* argv[]){
    std::ifstream fin(argv[1]);
    int f=(argc>2);
    std::ofstream fout("kokkos_realization/tmp.txt");
    Kokkos::initialize(argc, argv);{
        ll ord[256];
        ord['A']=1LL;
        ord['C']=2LL;
        ord['G']=3LL;
        ord['T']=4LL;
        std::string data;
        fin>>data;
        ll size=data.size();
        int len=(data.size()>10000?1000:data.size()/10);
        len=(len<3?3:len);
        std::vector<ll> powmod=pre(size);
        Kokkos::View<double*> freq("array", size);
        ll def_val=0;
        for(int i=0;i<len;++i) def_val=(def_val+ord[data[i]]*powmod[i])%mod;
        Kokkos::Timer timer;
        timer.reset();
        double st = timer.seconds();
        Kokkos::parallel_for( "yAx", size-len+1, KOKKOS_LAMBDA (int i) {
            ll res=0;
            ll hash=0;
            ll tmp_h=0;
#pragma omp simd
            for(int j=0;j<len;++j){
                hash=(hash+ord[data[i+j]]*powmod[j])%mod;
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
            for(int j=1;j<=size-len;++j){
                tmp_h=(tmp_h-(ord[data[j-1]]*powmod[j-1])%mod+(ord[data[len+j-1]]*powmod[len+j-1])%mod+mod)%mod;
                if (tmp_h==(hash*powmod[j])%mod){
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
        if (f) fout<<timer.seconds()-st<<" ";
    }
    Kokkos::finalize();
    return 0;
}
