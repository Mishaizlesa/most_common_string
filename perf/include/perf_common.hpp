#ifndef PERF_COMMON_HPP
#define PERF_COMMON_HPP


#include "algorithms.h"

#include <iostream>
#include <limits>
#include <random>
#include <vector>
#include <string>
#include <stdint.h>

#ifdef __riscv
static inline uint64_t __attribute__((__gnu_inline__, __always_inline__, __artificial__)) rdcycle(void)
{
    uint64_t dst;
    asm volatile ("csrrs %0, 0xc00, x0" : "=r" (dst));
    return dst;
}
#endif // __riscv

#ifdef __x86_64__
static __inline __attribute__((__gnu_inline__, __always_inline__, __artificial__)) uint64_t rdcycle(void){
#ifdef C_MSVC
  return __rdtsc();
#else
  uint64_t a, d;

  __asm__ __volatile__ ("rdtsc" : "=a" (a), "=d" (d));

  return ((uint64_t)a + ((uint64_t)d << 32));
#endif
}
#endif // __x86_64__

#endif // PERF_COMMON_H
