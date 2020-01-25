cmake_minimum_required(VERSION 3.4)
project(pwhd)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -mavx2")


if(${CMAKE_BUILD_TYPE} MATCHES "Release")
    message("Flags: ${CMAKE_CXX_FLAGS_RELEASE}")
endif()

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(BENCH_FILES xp-params.h bench/cli.cpp bench/cli.h bench/params.h bench/bench.cpp bench/bench.h)

set(HELPER_FILES
        utils/byteswap.h utils/helper.cpp utils/helper.h utils/testdata.h)

set(PWHD_FILES
        ${HELPER_FILES} ${BENCH_FILES}
        PwHammDistAlgorithmBase.h PwHammDistAlgorithmFactory.cpp PwHammDistAlgorithmFactory.h
        QuantizationPwHammDistAlgorithm.h
        )

add_executable(bench-binary bench/main-binary.cpp ${PWHD_FILES})
add_executable(bench-dataset bench/main-dataset.cpp ${PWHD_FILES} )

add_executable(bench-binary-win bench/main-binary.cpp ${PWHD_FILES})
add_executable(bench-dataset-win bench/main-dataset.cpp ${PWHD_FILES})
if((CYGWIN OR MINGW) AND ${CMAKE_BUILD_TYPE} MATCHES "Release")
    message("Windows static build: ${CMAKE_SYSTEM} ${CMAKE_SYSTEM_NAME}")
    target_link_libraries(bench-binary-win -static -static-libgcc -static-libstdc++)
    target_link_libraries(bench-dataset-win -static -static-libgcc -static-libstdc++)
endif()