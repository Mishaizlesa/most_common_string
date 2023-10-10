
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cuda.h>
#include <stdio.h>
#include <unordered_map>
typedef unsigned long long ll;


__global__ void search_kernel(ll* frequency, const short* data ,const ll len, const ll size) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i > size) return;
    int shift[64];
    for (int j = 0; j < 64; ++j) shift[j] = len - 2;
    int res = 0;
    int sh1;
    ll hash = 0;
    for (int j = 2; j <= len - 1; ++j) {
        int ind = (data[i + j - 2]) * 16 + (data[i + j - 1]) * 4 + (data[i + j]);
        if (j == len - 1) sh1 = shift[ind];
        shift[ind] = len - 1 - j;
    }

    if (!sh1) sh1 = 1;
    int j = len - 1;
    for (;;) {
        int sh = 1;
        while (sh && j <= size - len) {
            int ind = (data[j - 2]) * 16 + (data[j - 1]) * 4 + (data[j]);
            sh = shift[ind];
            j += sh;
        }
        if (j <= size - len) {
            int is_eq = 1;
            for (int k = 0; k < len; ++k) {
                if (data[i + k] != data[j - len + 1 + k]) {
                    is_eq = 0;
                    break;
                }
            }
            res += is_eq;
            j += sh1;
        }
        else {
            break;
        }
    }
    frequency[i] = res;
}

int main(int argc, char* argv[]) {
    std::unordered_map<char, char> mapSymbToCode = { {'A', (char)0}, {'C', (char)1}, {'G', (char)2}, {'T', (char)3} };
    float exe_milliseconds = 0;
    float copy_milliseconds = 0;
    cudaEvent_t start, exe_stop, copy_stop, start_copy;
    cudaEventCreate(&start);
    cudaEventCreate(&start_copy);
    cudaEventCreate(&exe_stop);
    cudaEventCreate(&copy_stop);

    ll len = std::atoi(argv[2]);
    std::ofstream fout("tmp.txt");
    FILE* data_file = fopen(argv[1], "rb");

    fseek(data_file, 0, SEEK_END);
    ll fsize = ftell(data_file);
    fseek(data_file, 0, SEEK_SET);
    char* data_h = (char*)malloc(fsize);


    fread(data_h, fsize, 1, data_file);

    short* enc_data_h = (short*)malloc(sizeof(short) * fsize);

    for (int i = 0; i < fsize; ++i) enc_data_h[i] = mapSymbToCode[data_h[i]];

    short* enc_data_d;
    cudaMalloc((void**)&enc_data_d, sizeof(short) * fsize);
    cudaMemcpy(enc_data_d, enc_data_h, fsize * sizeof(short), cudaMemcpyHostToDevice);
    ll* frequency_dev;
    cudaMalloc((void**)&frequency_dev, fsize*sizeof(ll) );

    ll* frequency_host = (ll*)malloc(fsize * sizeof(ll));
    cudaEventRecord(start,0);
    cudaEventRecord(start_copy, 0);
    search_kernel<< <(fsize - len + 256) / 256, 256 >> > (frequency_dev, enc_data_d, len, fsize);
    cudaDeviceSynchronize();

    cudaEventRecord(exe_stop,0);


    cudaEventSynchronize(exe_stop);

    cudaEventElapsedTime(&exe_milliseconds, start, exe_stop);
   
    cudaMemcpy(frequency_host, frequency_dev, fsize * sizeof(ll), cudaMemcpyDeviceToHost);


    cudaEventRecord(copy_stop, 0);
    cudaEventSynchronize(copy_stop);
    cudaEventElapsedTime(&copy_milliseconds, start_copy, copy_stop);
    std::cout<<exe_milliseconds<<" "<<copy_milliseconds;
}