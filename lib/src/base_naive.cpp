#include "algorithms.h"
#include "common_defs.h"
//#include "riscv_vector.h"
typedef long long ll;
extern void base_naive(std::vector<uint32_t>& freq ,const std::string& input_file, const uint32_t len_){
    std::ifstream fin(input_file);
    std::string data_;
    fin>>data_;
    ll size=data_.size();
    int len=len_;


    std::unordered_map<char,char>symbols_code{{'A',char(0)}, {'C',char(1)},{'G',char(2)},{'T',char(3)}};
    freq.resize(size);
    char data[size];
    for(int i=0;i<size;++i){
        data[i]=symbols_code[data_[i]];
    }
    int res;
    for (int i=0;i<size-len+1;++i){
        res=0;
        for(int j=0;j<size-len+1;++j){
            int is_eq=1;
            for(int k=0;k<len && is_eq;++k){
                if (data[i+k]!=data[j+k]) is_eq=0;
            }
            res+=is_eq;
        }
        freq[i] = res;
    }
    return;
}
