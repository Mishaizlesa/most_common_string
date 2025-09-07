cmake_minimum_required(VERSION 3.10)
project(most_common_string LANGUAGES CXX Fortran) # Добавляем Fortran

# Найти Intel SYCL
find_package(IntelSYCL REQUIRED)
# Базовая конфигурация для всех файлов
add_library(BaseConfiguration INTERFACE)
target_compile_options(BaseConfiguration
    INTERFACE
        -Wall
        -Wextra
)

if(BUILD_STATIC)
    target_link_options(BaseConfiguration INTERFACE -dynamic)
endif()

if(USE_OPENMP)
    target_compile_options(BaseConfiguration INTERFACE -qopenmp)
endif()

# Конфигурация архитектуры
if("${TARGET_ARCH}" STREQUAL "RISCV_GENERIC")
    target_compile_options(BaseConfiguration INTERFACE -march=rv64gc)
    target_compile_definitions(BaseConfiguration INTERFACE -DRISCV_GENERIC)
elseif("${TARGET_ARCH}" STREQUAL "RISCV_VECTOR")
    target_compile_options(BaseConfiguration INTERFACE -march=rv64gcv)
    target_compile_definitions(BaseConfiguration INTERFACE -DRISCV_VECTOR)
elseif("${TARGET_ARCH}" STREQUAL "X86")
    target_compile_definitions(BaseConfiguration INTERFACE -DX86)
    ##target_compile_options(BaseConfiguration INTERFACE -mavx512bw)
    ##target_compile_options(BaseConfiguration INTERFACE -mavx512f)
    ##target_compile_options(BaseConfiguration INTERFACE -mavx512vl)
    ##target_compile_options(BaseConfiguration INTERFACE -mavx512dq)
    ##target_compile_options(BaseConfiguration INTERFACE -mavx2)
##
    ##target_compile_options(BaseConfiguration INTERFACE -fsycl)
    ##target_compile_options(BaseConfiguration INTERFACE -std=c++17)
else()
    message(FATAL_ERROR "Unsupported TARGET_ARCH")
endif()

# Добавляем SYCL и C++17
target_compile_options(BaseConfiguration INTERFACE -fsycl -std=c++17 -m64)

# Конфигурация оптимизации
add_library(BaseOptConfiguration INTERFACE)
if(BUILD_TYPE STREQUAL "Release")
    target_compile_options(BaseOptConfiguration INTERFACE -O3)
else()
    target_compile_options(BaseOptConfiguration
        INTERFACE
            -O0
            -g
            -DNDEBUG
    )
endif()

# Общая конфигурация
add_library(CommonConfiguration INTERFACE)
target_link_libraries(CommonConfiguration INTERFACE BaseConfiguration BaseOptConfiguration)