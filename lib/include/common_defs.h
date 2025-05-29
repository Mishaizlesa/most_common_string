#ifndef COMMON_H
#define COMMON_H

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <omp.h>
#include <cstdio>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <unordered_map>
//#include "riscv_vector.h"
//typedef __rvv_bool64_t vbool64_t;
//typedef __rvv_bool32_t vbool32_t;
//typedef __rvv_bool16_t vbool16_t;
//typedef __rvv_bool8_t vbool8_t;
//typedef __rvv_bool4_t vbool4_t;
//typedef __rvv_bool2_t vbool2_t;
//typedef __rvv_bool1_t vbool1_t;
//typedef __rvv_int8mf8_t vint8mf8_t;
//typedef __rvv_uint8mf8_t vuint8mf8_t;
//typedef __rvv_int8mf4_t vint8mf4_t;
//typedef __rvv_uint8mf4_t vuint8mf4_t;
//typedef __rvv_int8mf2_t vint8mf2_t;
//typedef __rvv_uint8mf2_t vuint8mf2_t;
//typedef __rvv_int8m1_t vint8m1_t;
//typedef __rvv_uint8m1_t vuint8m1_t;
//typedef __rvv_int8m2_t vint8m2_t;
//typedef __rvv_uint8m2_t vuint8m2_t;
//typedef __rvv_int8m4_t vint8m4_t;
//typedef __rvv_uint8m4_t vuint8m4_t;
//typedef __rvv_int8m8_t vint8m8_t;
//typedef __rvv_uint8m8_t vuint8m8_t;
//typedef __rvv_int16mf4_t vint16mf4_t;
//typedef __rvv_uint16mf4_t vuint16mf4_t;
//typedef __rvv_int16mf2_t vint16mf2_t;
//typedef __rvv_uint16mf2_t vuint16mf2_t;
//typedef __rvv_int16m1_t vint16m1_t;
//typedef __rvv_uint16m1_t vuint16m1_t;
//typedef __rvv_int16m2_t vint16m2_t;
//typedef __rvv_uint16m2_t vuint16m2_t;
//typedef __rvv_int16m4_t vint16m4_t;
//typedef __rvv_uint16m4_t vuint16m4_t;
//typedef __rvv_int16m8_t vint16m8_t;
//typedef __rvv_uint16m8_t vuint16m8_t;
//typedef __rvv_int32mf2_t vint32mf2_t;
//typedef __rvv_uint32mf2_t vuint32mf2_t;
//typedef __rvv_int32m1_t vint32m1_t;
//typedef __rvv_uint32m1_t vuint32m1_t;
//typedef __rvv_int32m2_t vint32m2_t;
//typedef __rvv_uint32m2_t vuint32m2_t;
//typedef __rvv_int32m4_t vint32m4_t;
//typedef __rvv_uint32m4_t vuint32m4_t;
//typedef __rvv_int32m8_t vint32m8_t;
//typedef __rvv_uint32m8_t vuint32m8_t;
//typedef __rvv_int64m1_t vint64m1_t;
//typedef __rvv_uint64m1_t vuint64m1_t;
//typedef __rvv_int64m2_t vint64m2_t;
//typedef __rvv_uint64m2_t vuint64m2_t;
//typedef __rvv_int64m4_t vint64m4_t;
//typedef __rvv_uint64m4_t vuint64m4_t;
//typedef __rvv_int64m8_t vint64m8_t;
//typedef __rvv_uint64m8_t vuint64m8_t;
#endif // COMMON_H