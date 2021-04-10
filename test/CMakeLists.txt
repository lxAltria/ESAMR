set (EVA_INCLUDES "../external/compression_utils/include")

add_executable (test_encoder test_encoder.cpp)
target_include_directories(test_encoder PRIVATE ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES})
target_link_libraries(test_encoder ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB})

add_executable (test_decomposer test_decomposer.cpp)
target_include_directories(test_decomposer PRIVATE ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES})
target_link_libraries(test_decomposer ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB})

add_executable (test_interleaver test_interleaver.cpp)
target_include_directories(test_interleaver PRIVATE ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES})
target_link_libraries(test_interleaver ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB})

add_executable (test_error_collector test_error_collector.cpp)
target_include_directories(test_error_collector PRIVATE ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES})
target_link_libraries(test_error_collector ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB})

add_executable (test_refactor test_refactor.cpp)
target_include_directories(test_refactor PRIVATE ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES} ${MPI_INCLUDE_PATH})
target_link_libraries(test_refactor ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB} ${MPI_LIBRARIES} adios2::adios2)

add_executable (test_reconstructor test_reconstructor.cpp)
target_include_directories(test_reconstructor PRIVATE ${EVA_INCLUDES} ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES} ${MPI_INCLUDE_PATH})
target_link_libraries(test_reconstructor ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB} ${MPI_LIBRARIES} adios2::adios2)
