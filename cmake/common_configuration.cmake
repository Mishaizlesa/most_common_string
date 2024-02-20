
# Base configuration for all files in project.
# (including custom configurations in lib/CMakeLists.txt)
add_library(BaseConfiguration INTERFACE)

target_compile_options(BaseConfiguration
INTERFACE
    -Wall
    -Wextra
)

if(BUILD_STATIC)
    target_link_options(BaseConfiguration INTERFACE -static)
endif()

if(USE_OPENMP)
    target_compile_options(BaseConfiguration INTERFACE -fopenmp)
    target_link_options(BaseConfiguration INTERFACE -fopenmp)
    target_compile_definitions(BaseConfiguration INTERFACE -DUSE_OPENMP)
endif()

if("${TARGET_ARCH}" STREQUAL "RISCV_GENERIC")
    target_compile_options(BaseConfiguration INTERFACE -march=rv64gc)
    target_compile_definitions(BaseConfiguration INTERFACE -DRISCV_GENERIC)
elseif("${TARGET_ARCH}" STREQUAL "RISCV_VECTOR")
    target_compile_options(BaseConfiguration INTERFACE -march=rv64gcv0p7)
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
