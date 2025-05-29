
# Base configuration for all files in project.
# (including custom configurations in lib/CMakeLists.txt)
add_library(BaseConfiguration INTERFACE)

target_compile_options(BaseConfiguration
INTERFACE
    -Wall
    -Wextra
)

set(KOKKOS_ROOT "/home/mixa/kokkos-avx512" CACHE PATH "Path to Kokkos installation")

# Поиск Kokkos
find_package(Kokkos REQUIRED 
    PATHS ${KOKKOS_ROOT}
    PATH_SUFFIXES lib/cmake/Kokkos lib64/cmake/Kokkos
    NO_DEFAULT_PATH
)

target_compile_options(BaseConfiguration INTERFACE ${RISCV_FLAGS})
target_include_directories(BaseConfiguration INTERFACE 
    ${Kokkos_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}/lib/include
    ${PROJECT_SOURCE_DIR}/perf/include
)
target_link_libraries(BaseConfiguration INTERFACE Kokkos::kokkos)


if(BUILD_STATIC)
    target_link_options(BaseConfiguration INTERFACE -static)
endif()

if(USE_OPENMP)
   # target_compile_options(BaseConfiguration INTERFACE -fopenmp)
    # Явно указываем статическую библиотеку
    target_link_libraries(BaseConfiguration INTERFACE 
        /home/mixa/newtoolkit/sc-devtoolkit/riscv-gcc/riscv64-unknown-linux-gnu/lib64/lp64d/libgomp.a
    )
endif()

if("${TARGET_ARCH}" STREQUAL "RISCV_GENERIC")
    target_compile_options(BaseConfiguration INTERFACE -march=rv64gc)
    target_compile_definitions(BaseConfiguration INTERFACE -DRISCV_GENERIC)
elseif("${TARGET_ARCH}" STREQUAL "RISCV_VECTOR")
    target_compile_options(BaseConfiguration INTERFACE -march=skylake-avx512)
    target_compile_definitions(BaseConfiguration INTERFACE -DRISCV_VECTOR)
elseif("${TARGET_ARCH}" STREQUAL "X86")
target_compile_definitions(BaseConfiguration INTERFACE -DX86)
else()
    message(ERROR "Unsupported TARGET_ARCH")
endif()

# Base configuration with optimization flags
# Used for all files except custom configuration in lib/CMakeLists.txt and benchmarks
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

# Common configuration.
# Used for all files except custom configuration in lib/CMakeLists.txt and benchmarks
add_library(CommonConfiguration INTERFACE)
target_link_libraries(CommonConfiguration INTERFACE BaseConfiguration BaseOptConfiguration)
