#include <Kokkos_Core.hpp>
#include <Kokkos_Atomic.hpp>
#include <omp.h>
#include <cstdio>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
typedef long long ll;
const ll p=31;
int main(int argc, char* argv[]){
    std::ifstream fin(argv[1]);
    int f=argv[2][0]-'0';
    std::ofstream fout("tmp.txt");
    Kokkos::initialize(argc, argv);{
        Kokkos::Timer timer;
        std::string data_;
        fin>>data_;
        ll size=data_.size();
        int len=(data_.size()>10000?1000:data_.size()/10);
        len=(len<3?3:len);
        Kokkos::View<double*,Kokkos::SharedSpace> freq("array", size);
        Kokkos::View<ll*,Kokkos::SharedSpace> dp_h("dp_h", size+1);
        Kokkos::View<char*,Kokkos::SharedSpace> data("device_string",size);
        for(int i=0;i<size;++i){
            data[i]=data_[i];
        }
        //std::vector<ll>dp_h(size+1);
        dp_h[0]=0;
        timer.reset();
        double st=timer.seconds();
        ll mod=1e9+7;
        ll powmod=1;
        for(int i=0;i<size;++i){
            if (data_[i]=='A') dp_h[i+1]=(dp_h[i]+powmod)%mod;
            else if (data_[i]=='C') dp_h[i+1]=(dp_h[i]+2LL*powmod)%mod;
            else if (data_[i]=='G') dp_h[i+1]=(dp_h[i]+3LL*powmod)%mod;
            else  dp_h[i+1]=(dp_h[i]+4LL*powmod)%mod;
            powmod*=p;
            powmod%=mod;
        }
        Kokkos::parallel_for( "yAx", size-len+1, KOKKOS_LAMBDA (int i) {
            int res=0;
            ll hash=0;
            ll mod=1e9+7;
            ll p=31;
            ll powmod=1;
            for(int j=0;j<len;++j){
                if (data[i+j]=='A') hash=(hash+powmod)%mod;
                else if (data[i+j]=='C') hash=(hash+2LL*powmod)%mod;
                else if (data[i+j]=='G') hash=(hash+3LL*powmod)%mod;
                else  hash=(hash+4LL*powmod)%mod;
                powmod*=p;
                powmod%=mod;
            }
            powmod=1;
            for(int j=0;j<=size-len;++j){
                ll val_h=(dp_h[j+len]+mod-dp_h[j])%mod;
                if (val_h==(hash*powmod)%mod){
                    int is_eq=1;
                    for(int k=0;k<len;++k){
                        if (data[i+k]!=data[j+k]){
                            is_eq=0;break;
                        }
                    }
                    res+=is_eq;
                }
                powmod*=p;
                powmod%=mod;
            }
            freq(i)=res;
        });
        if (f==1){
            fout<<timer.seconds()-st<<" ";
        }else if (f==2){
            for(int i=0;i<=size-len;++i){
                fout<<freq[i];
            }
        }
    }
    Kokkos::finalize();
    return 0;
}
