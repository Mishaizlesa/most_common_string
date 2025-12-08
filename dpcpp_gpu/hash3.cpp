#include <CL/sycl.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <chrono>

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

    char* data = sycl::malloc_shared<char>(size, q);
    int* freq = sycl::malloc_shared<int>(size, q);

    // Initialize data from file
    for (int i = 0; i < size; ++i) {
        data[i] = data_[i];
    }

    auto start = std::chrono::high_resolution_clock::now();
    q.submit([&](sycl::handler& h) {
        h.parallel_for(sycl::range<1>(size - len + 1), [=](sycl::id<1> i) {
            int shift[400];
            for (int j = 0; j < 400; ++j) shift[j] = len - 2;
            int res = 0;
            int sh1;
            ll hash = 0;
            
            for (int j = 2; j <= len - 1; ++j) {
                int ind = (data[i + j - 2] - 'A') * 16 + (data[i + j - 1] - 'A') * 4 + (data[i + j] - 'A');
                if (j == len - 1) sh1 = shift[ind];
                shift[ind] = len - 1 - j;
            }
            
            if (!sh1) sh1 = 1;
            int j = len - 1;
            
            for (;;) {
                int sh = 1;
                while (sh && j < size) {
                    int ind = (data[j - 2] - 'A') * 16 + (data[j - 1] - 'A') * 4 + (data[j] - 'A');
                    sh = shift[ind];
                    j += sh;
                }
                if (j < size) {
                    int is_eq = 1;
                    for (int k = 0; k < len; ++k) {
                        if (data[i + k] != data[j - len + 1 + k]) {
                            is_eq = 0;
                            break;
                        }
                    }
                    res += is_eq;
                    j += sh1;
                } else {
                    break;
                }
            }
            freq[i] = res;
        });
    }).wait();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    if (f == 1) { 
        std::cout << elapsed.count() << " ";
    }
    if (f == 2) {
        for (int i = 0; i <= size - len; ++i) {
            fout << freq[i];
        }
    } else if (f == 3) {
        int mmax = 0;
        for (int i = 0; i < size; ++i) {
            mmax = (mmax > freq[i] ? mmax : freq[i]);
        }
        fout << mmax;
    }

    sycl::free(data, q);
    sycl::free(freq, q);

    return 0;
}