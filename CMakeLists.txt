cmake_minimum_required(VERSION 3.16)
project(tt_a3_operator_credit CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# ISO C++ (-std=c++20, not gnu++20) and no fused multiply-add contraction: an FMA rounds once where a
# multiply and an add round twice, so a compiler that fuses (GCC does by default on arm64 in GNU
# mode) would change the floating-point trajectories and the records would not reproduce across
# architectures.  On x86-64 without -mfma nothing is fused either way.
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffp-contract=off")
find_package(Threads REQUIRED)
file(GLOB CORE_SRC ${CMAKE_SOURCE_DIR}/src/core/*.cpp)
add_executable(experiments src/main.cpp ${CORE_SRC})
target_include_directories(experiments PRIVATE src)
target_link_libraries(experiments PRIVATE Threads::Threads)
