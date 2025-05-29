

# Project build settings
ENABLE_TEST="ON" # ON or OFF
ENABLE_BENCH="ON" # ON or OFF
USE_OPENMP="ON" # ON or OFF
TARGET_ARCH="RISCV_VECTOR" # RISCV_GENERIC or RISCV_VECTOR
BUILD_TYPE="Release" # Release or Debug
BUILD_FOLDER="_build_x86"
# Path to C and C++ compilers
C_COMPILER_PATH="icx"
CXX_COMPILER_PATH="icx"

# Clear build folder
# Comment if not needed
rm -rf $BUILD_FOLDER

# Confuigure project
cmake CMakeLists.txt \
 -DENABLE_TEST=$ENABLE_TEST \
 -DENABLE_BENCH=$ENABLE_BENCH \
 -DTARGET_ARCH=$TARGET_ARCH \
 -B"$BUILD_FOLDER" \
 -DCMAKE_C_COMPILER=$C_COMPILER_PATH \
 -DCMAKE_CXX_COMPILER=$CXX_COMPILER_PATH \
 -DBUILD_TYPE=$BUILD_TYPE \
 -DBUILD_STATIC=ON

# Build project
cmake --build  _build_x86