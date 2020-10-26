cmake_minimum_required (VERSION 3.0)
project(MDR VERSION 0.0.0.1)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_FLAGS "-g -O3")

set (MGARDx_INCLUDES "/Users/xin/github/test/MGARDx/install/include")

find_library(ZSTD_LIB zstd HINTS "/Users/xin/github/SZ_master/build/zstd")
find_library(SZ3_LIB SZ3 HINTS "/Users/xin/github/SZ_refactory/build/src")
set (ZSTD_INCLUDES "/Users/xin/github/SZ_master/zstd")
set (SZ3_INCLUDES "/Users/xin/github/SZ_refactory/include")

add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(${PROJECT_NAME} INTERFACE include)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include)
add_subdirectory (test)
