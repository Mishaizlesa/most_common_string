#include <CL/sycl.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <chrono>
#include <vector>
#include <algorithm>
#include <unordered_map>

typedef long long ll;

int main(int argc, char* argv[]) {
    std::ifstream fin(argv[1]);
    int f = argv[2][0] - '0';
    std::ofstream fout("tmp.txt");
    
    std::string data_;
    fin >> data_;
    ll size = data_.size();
    int len = std::atoi(argv[3]);
    len = (len < 3 ? 3 : len);

    sycl::queue q(sycl::gpu_selector_v);
    std::cout << "Running on device: " << q.get_device().get_info<sycl::info::device::name>() << "\n";

    std::unordered_map<char, char> mapSymbToCode = {
        {'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}
    };
    
    std::vector<char> coded_data(size);
    for (ll i = 0; i < size; ++i) {
        coded_data[i] = mapSymbToCode[data_[i]];
    }

    char* data_device = sycl::malloc_shared<char>(size, q);
    int* freq_device = sycl::malloc_shared<int>(size, q);
    
    q.memcpy(data_device, coded_data.data(), size * sizeof(char)).wait();
    q.memset(freq_device, 0, size * sizeof(int)).wait();

    size_t global_size = size - len + 1;
    size_t local_size = 256;

    auto start = std::chrono::high_resolution_clock::now();
    
    q.submit([&](sycl::handler& h) {
        h.parallel_for(sycl::nd_range<1>(
            sycl::range<1>((global_size + local_size - 1) / local_size * local_size),
            sycl::range<1>(local_size)
        ), [=](sycl::nd_item<1> item) {
            size_t i = item.get_global_id(0);
            if (i >= global_size) return;

            uint64_t hash_pattern = 0;
            uint64_t hash_text = 0;
            uint64_t p = 2;
            uint64_t powmod = 1;
            int res = 0;
            
            for (int j = 0; j < len; ++j) {
                hash_pattern = hash_pattern * p + data_device[i + j];
                hash_text = hash_text * p + data_device[j];
                powmod *= p;
            }
            
            for (size_t j = 0; j < global_size; ++j) {
                if (hash_text == hash_pattern) {
                    int is_eq = 1;
                    for (int k = 0; k < len; ++k) {
                        if (data_device[i + k] != data_device[j + k]) {
                            is_eq = 0;
                            break;
                        }
                    }
                    res += is_eq;
                }
                if (j < global_size - 1) {
                    hash_text = hash_text * p - data_device[j] * powmod + data_device[j + len];
                }
            }
            
            freq_device[i] = res;
        });
    }).wait();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::vector<int> freq_host(size);
    q.memcpy(freq_host.data(), freq_device, size * sizeof(int)).wait();

    if (f == 1) { 
        std::cout << elapsed.count() << " ";
    }
    if (f == 2) {
        for (int i = 0; i <= size - len; ++i) {
            fout << freq_host[i];
        }
    } else if (f == 3) {
        int mmax = 0;
        for (int i = 0; i <= size - len; ++i) {
            mmax = (mmax > freq_host[i] ? mmax : freq_host[i]);
        }
        fout << mmax;
    }

    sycl::free(data_device, q);
    sycl::free(freq_device, q);

    return 0;
}