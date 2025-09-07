#include <sycl/sycl.hpp>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <unordered_map>
#include "common_defs.h"

extern "C" void hash3_scalar(std::vector<uint32_t>& freq, const std::string& input_file, 
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
      int32_t shift[64];
      for (int k = 0; k < 64; ++k) {
        shift[k] = len_ - 2;
      }
      for (int j = 2; j < len_ - 1; ++j) {
        int ind = data_acc[i + j - 2] * 16 + data_acc[i + j - 1] * 4 + data_acc[i + j];
        shift[ind] = len_ - 1 - j;
      }
      int ind = data_acc[i + len_ - 3] * 16 + data_acc[i + len_ - 2] * 4 + data_acc[i + len_ - 1];
      int sh1 = shift[ind];
      shift[ind] = 0;
      if (sh1 == 0) sh1 = 1;
      int j = len_ - 1;
      while (true) {
        int sh = 1;
        while (sh != 0 && j < N) {
          ind = data_acc[j - 2] * 16 + data_acc[j - 1] * 4 + data_acc[j];
          sh = shift[ind];
          j += sh;
        }
        if (j < N) {
          bool is_eq = true;
          for (uint32_t k = 0; k < len_; ++k) {
            if (data_acc[i + k] != data_acc[j - len_ + 1 + k]) {
              is_eq = false;
              break;
            }
          }
          if (is_eq) ++res;
          j += sh1;
        } else {
          break;
        }
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