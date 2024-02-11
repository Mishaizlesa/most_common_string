#include <omp.h>
#include <iostream>
#include <string>
#include <unordered_map>
#include <fstream>
typedef long long ll;
int main(int argc, char* argv[]){
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



    #pragma omp parallel for shared(data, freq)
        for (int i=0;i<size-len+1;++i){
            for(int j=0;j<size-len+1;++j){
                int is_eq=1;
                for(int k=0;k<len && is_eq;++k){
                    if (data[i+k]!=data[j+k]) is_eq=0;
                }
                freq[i]+=is_eq;
            }
        }
    
    


    if (f==1){
        //fout<<timer.seconds()<<" ";
    }
    if (f==2){
        std::ofstream fout(argv[3]);
        for(int i=0;i<size-len+1;++i){
            fout<<freq[i];
        }
    }
    return 0;
}
