add_executable (test_refactor test_refactor.cpp)
target_include_directories(test_refactor PRIVATE ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES})
target_link_libraries(test_refactor ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB})

add_executable (test_retrieval test_retrieval.cpp)
target_include_directories(test_retrieval PRIVATE ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES})
target_link_libraries(test_retrieval ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB})

# find_library(ZFP_LIB zfp HINTS "/Users/xin/github/zfp/build/lib")
# set (ZFP_INC "/Users/xin/github/zfp/include")

# add_executable (test_bitplane_performance test_bitplane_performance.cpp)
# target_include_directories(test_bitplane_performance PRIVATE ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES} ${ZFP_INC})
# target_link_libraries(test_bitplane_performance ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB} ${ZFP_LIB})
