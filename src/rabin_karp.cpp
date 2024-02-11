#include <omp.h>
#include <cstdio>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <unordered_map>
typedef long long ll;

int main(int argc, char* argv[]) {
    std::ifstream fin(argv[1]);
    int f = std::stoi(argv[2]);
    int len = std::stoi(argv[3]);
    std::ofstream fout("tmp.txt");

    std::string data_;
    fin >> data_;
    ll size = data_.size();

    //std::cout << size << " " << len << std::endl;
    int freq[size];
    char data[size+1];
    data[size] = 0; 
    std::unordered_map<char, char> mapSymbToCode = { {'A', (char)0}, {'C', (char)1}, {'G', (char)2}, {'T', (char)3} };
    // prepare data
    for (int i = 0; i < size; ++i) {
        data[i] = mapSymbToCode[data_[i]];
    }


    //std::cout << timer.seconds() - st << std::endl;
    #pragma omp parallel for shared(data, freq)
    for (int i=0;i<size-len+1;++i){
        //for (int i = 0; i < size - len + 1; i++) {
        int res = 0;
        ll hash_pattern = 0;
        ll hash_text = 0;
        ll mod = 1LL<<31 - 1;
        ll p = 2;
        ll powmod = 1;
        for (int j = 0; j < len; ++j) {
            hash_pattern = hash_pattern * p + data[i + j];
            hash_text = hash_text * p + data[j];
            powmod = powmod * p;
    powmod%=mod;
    hash_pattern %= mod;
            hash_text %= mod;
        }
        for (int j = 0; j < size - len + 1; ++j) {
            if (hash_text == hash_pattern) {
                int is_eq = 1;
                for (int k = 0; k < len && is_eq; ++k) {
                    if (data[i + k] != data[j + k]) {
            is_eq=0;
                    }
                }
                res += is_eq;
            }
            hash_text = (hash_text * p - data[j] * powmod + data[j + len]) % mod;
        }
        freq[i] = res;
    }
    if (f == 1) {
        fout << 0 << " ";
    }
    else if (f == 2) {
        for (int i = 0; i <= size - len; ++i) {
            fout << freq[i];
        }
    }
    return 0;
}