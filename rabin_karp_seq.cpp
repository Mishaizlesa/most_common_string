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
    int f=(argc>2);
    std::ofstream fout("tmp.txt");
    Kokkos::Timer timer;
    std::string data;
    fin>>data;
    ll size=data.size();
    int len=(data.size()>10000?1000:data.size()/10);
    len=(len<3?3:len);
    std::vector<ll>freq(size);
    std::vector<ll>dp_h(size+1);
    dp_h[0]=0;
    timer.reset();
    double st=timer.seconds();
    ll powmod=1;
    for(int i=0;i<size;++i){
        if (data[i]=='A') dp_h[i+1]=(dp_h[i]+powmod)%mod;
        else if (data[i]=='C') dp_h[i+1]=(dp_h[i]+2LL*powmod)%mod;
        else if (data[i]=='G') dp_h[i+1]=(dp_h[i]+3LL*powmod)%mod;
        else  dp_h[i+1]=(dp_h[i]+4LL*powmod)%mod;
        powmod*=p;
        powmod%=mod;
    }
    for(int i=0;i<=size-len;++i){
        int res=0;
        ll hash=0;
        powmod=1;
#pragma omp simd
        for(int j=0;j<len;++j){
            if (data[i+j]=='A') hash=(hash+powmod)%mod;
            else if (data[i+j]=='C') hash=(hash+2LL*powmod)%mod;
            else if (data[i+j]=='G') hash=(hash+3LL*powmod)%mod;
            else  hash=(hash+4LL*powmod)%mod;
            powmod*=31;
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
            powmod*=31;
            powmod%=mod;
        }
        freq[i]=res;
    }
    if (f) fout<<timer.seconds()-st<<" ";
    return 0;
}
