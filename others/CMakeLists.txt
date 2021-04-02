set(ZFP_INCLUDES "/gpfs/alpine/csc143/proj-shared/gongq/andes/MReduction/xin/zfp/include")
find_library(ZFP_LIB zfp HINTS "/gpfs/alpine/csc143/proj-shared/gongq/andes/MReduction/xin/zfp/lib")
set(SZ_INCLUDES "/gpfs/alpine/csc143/proj-shared/gongq/andes/MReduction/xin/SZ/install/include/sz")
find_library(SZ_LIB SZ HINTS "/gpfs/alpine/csc143/proj-shared/gongq/andes/MReduction/xin/SZ/install/lib64")

add_executable (parallel_sz parallel_sz.cpp)
target_include_directories(parallel_sz PRIVATE ${SZ_INCLUDES} ${ZSTD_INCLUDES})
target_link_libraries(parallel_sz ${SZ_LIB} ${MPI_LIBRARIES} adios2::adios2)

add_executable (parallel_zfp parallel_zfp.cpp)
target_include_directories(parallel_zfp PRIVATE ${ZFP_INCLUDES})
target_link_libraries(parallel_zfp ${ZFP_LIB} ${MPI_LIBRARIES} adios2::adios2)

find_library(SZ_CPP_LIB sz_cpp HINTS "/gpfs/alpine/csc143/proj-shared/gongq/andes/MReduction/xin/Multiprecision-data-refactoring/external/meta_compressor/build/src")
set (SZ_CPP_INC "/gpfs/alpine/csc143/proj-shared/gongq/andes/MReduction/xin/Multiprecision-data-refactoring/external/meta_compressor/build/include")
add_executable (parallel_mgard parallel_mgard.cpp)
target_include_directories(parallel_mgard PRIVATE ${SZ3_INCLUDES} ${SZ_CPP_INC} ${ZSTD_INCLUDES} ${MGARDx_INCLUDES})
target_link_libraries(parallel_mgard ${ZSTD_LIB} ${SZ3_LIB} ${SZ_CPP_LIB} ${MPI_LIBRARIES} adios2::adios2)

