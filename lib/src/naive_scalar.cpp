#include <sycl/sycl.hpp>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <unordered_map>
#include "common_defs.h"

extern "C" void naive_scalar(std::vector<uint32_t>& freq, const std::string& input_file, 
                                 const uint32_t len_, const bool perf_collect) {
  std::ifstream fin(input_file);
  std::string data_str;
  fin >> data_str;
  size_t N = data_str.size();
  size_t M = N - len_ + 1;
  freq.resize(N, 0);

  std::vector<int8_t> data(N);
  std::unordered_map<char, int8_t> symbols_code = {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};
  for (size_t i = 0; i < N; ++i) {
    data[i] = symbols_code[data_str[i]];
  }

  auto start = std::chrono::high_resolution_clock::now();

  sycl::queue q{sycl::cpu_selector_v};

  sycl::buffer<int8_t> data_buf(data.data(), N);
  sycl::buffer<uint32_t> freq_buf(freq.data(), N);

  q.submit([&](sycl::handler& h) {
    auto data_acc = data_buf.get_access<sycl::access::mode::read>(h);
    auto freq_acc = freq_buf.get_access<sycl::access::mode::write>(h);

    h.parallel_for(M, [=](sycl::id<1> idx) {
      size_t i = idx[0];
      uint32_t res = 0;
      for (size_t j = 0; j < M; ++j) {
        bool is_eq = true;
        for (uint32_t k = 0; k < len_; ++k) {
          if (data_acc[i + k] != data_acc[j + k]) {
            is_eq = false;
            break;
          }
        }
        if (is_eq) ++res;
      }
      freq_acc[i] = res;
    });
  });
  q.wait();

  auto stop = std::chrono::high_resolution_clock::now();
  double elapsed = std::chrono::duration<double>(stop - start).count();

  if (perf_collect) {
    std::cout << N << " " << len_ << " " << elapsed << "\n";
  }
}