#include "perf_common.hpp"
int main() {
   //std::string path = "/home/mixa/most_common_string_kokkos/genome_samples/s103258.txt";
    std::string path ="genome_samples/s43794.txt";
    std::vector<uint32_t>freq;
    for (int i=4;i<=1024;i*=2){
  
        hash3(freq, path, i, true);
    }

    return 0;
}
