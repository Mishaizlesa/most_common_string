#include <CL/sycl.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <chrono>

typedef long long ll;

int main(int argc, char* argv[]) {

    std::ifstream fin(argv[1]);
    int f = argv[2][0] - '0';
    int len = std::atoi(argv[3]);
    
    std::string data_;
    fin >> data_;
    ll size = data_.size();
    ll valid_range = size - len + 1;

    sycl::queue q(sycl::gpu_selector_v);
    std::cout << q.get_device().get_info<sycl::info::device::name>() << "\n";


    char* data_device = sycl::malloc_device<char>(size, q);
    int* freq_device = sycl::malloc_device<int>(valid_range, q);
    
    q.memcpy(data_device, data_.data(), size * sizeof(char)).wait();
    q.memset(freq_device, 0, valid_range * sizeof(int)).wait();

    constexpr size_t local_size = 256;
    const size_t global_size = ((valid_range + local_size - 1) / local_size) * local_size;

    auto start = std::chrono::high_resolution_clock::now();

    q.submit([&](sycl::handler& h) {
        h.parallel_for(sycl::nd_range<1>(global_size, local_size), 
        [=](sycl::nd_item<1> item) {
            ll i = item.get_global_id(0);
            if (i >= valid_range) return;

            int count = 0;
            
            for (ll j = 0; j < valid_range; j++) {
                bool match = true;

                for (int k = 0; k < len; k += 4) {
                    if (k + 3 < len) {
                        if (data_device[i + k] != data_device[j + k] ||
                            data_device[i + k + 1] != data_device[j + k + 1] ||
                            data_device[i + k + 2] != data_device[j + k + 2] ||
                            data_device[i + k + 3] != data_device[j + k + 3]) {
                            match = false;
                            break;
                        }
                    } else {
                        for (int r = k; r < len; ++r) {
                            if (data_device[i + r] != data_device[j + r]) {
                                match = false;
                                break;
                            }
                        }
                        break;
                    }
                }
                count += match ? 1 : 0;
            }
            
            freq_device[i] = count;
        });
    }).wait();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::vector<int> freq_host(valid_range);
    q.memcpy(freq_host.data(), freq_device, valid_range * sizeof(int)).wait();

    if (f == 1) {
        std::cout << elapsed.count() << std::endl;
    } else if (f == 2) {
        std::ofstream fout(argv[3]);
        for (ll i = 0; i < valid_range; ++i) {
            fout << freq_host[i] << " ";
        }
    }

    sycl::free(data_device, q);
    sycl::free(freq_device, q);

    return 0;
}