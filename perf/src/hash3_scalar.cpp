#include "perf_common.hpp"
int main() {
   //std::string path = "/home/mixa/most_common_string_kokkos/genome_samples/s103258.txt";
    std::string path ="genome_samples/s201216.txt";
    std::vector<uint32_t>freq;
    for (int i=128;i<=128;i*=2){
  
        hash3_scalar(freq, path, i, true);
    }

    return 0;
}
