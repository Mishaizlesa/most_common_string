#include <omp.h>
#include <cstdio>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <unordered_map>
typedef long long ll;
int main(int argc, char* argv[]) {
     std::ifstream fin(argv[1]);
    int f=argv[2][0]-'0';
    std::ofstream fout("tmp.txt");
    std::string data_;
    fin>>data_;


    ll size=data_.size();
    int len=std::atoi(argv[3]);


    std::unordered_map<char,char>symbols_code{{'A',char(0)}, {'C',char(1)},{'G',char(2)},{'T',char(3)}};
    int freq[size];
    char data[size];
    for(int i=0;i<size;++i){
        data[i]=symbols_code[data_[i]];
    }



#pragma omp parallel for shared(freq)
    for(int i=0;i<=size-len;++i){
        int res=0;
        int sh1;
        std::vector<int>shift(64,len-2);
        ll hash=0;
        for(int j=2;j<=len-1;++j){
            int ind=data[i+j-2]*16+data[i+j-1]*4+data[i+j];
            if (j==len-1) sh1=shift[ind];
            shift[ind]=len-1-j;
        }
        
        if (!sh1) sh1=1;
        int j=len-1;
        for(;;){
            int sh=1;
            while (sh && j<size) {
                int ind=data[j-2]*16+data[j-1]*4+data[j];
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
    }
   // if (f==1) fout<<timer.seconds()-st<<" ";
    if(f==2){
        for(int i=0;i<=size-len;++i){
            fout<<freq[i];
        }
    }
    return 0;
}
