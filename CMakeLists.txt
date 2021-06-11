cmake_minimum_required (VERSION 3.0)
project(MDR VERSION 0.0.0.1)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_FLAGS "-g -O3")

set (MGARDx_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/external/MGARDx/install/include")

find_library(ZSTD_LIB zstd HINTS "${CMAKE_CURRENT_SOURCE_DIR}/external/SZ/install/lib")
find_library(SZ3_LIB SZ3 HINTS "${CMAKE_CURRENT_SOURCE_DIR}/external/SZ3/build/src")
find_library(MGARD_LIB mgard HINTS "${CMAKE_CURRENT_SOURCE_DIR}/external/MGARD/install/lib")
set (ZSTD_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/external/SZ/install/include")
set (SZ3_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/external/SZ3/include")
set (MGARD_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/external/MGARD/install/include")

add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(${PROJECT_NAME} INTERFACE include)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include)
add_subdirectory (test)
