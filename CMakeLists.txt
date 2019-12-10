cmake_minimum_required(VERSION 3.13)
project(dbpit)

set(CMAKE_CXX_STANDARD 14)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--stack,1000000000")
set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -mavx2")


if(${CMAKE_BUILD_TYPE} MATCHES "Release")
    message("Flags: ${CMAKE_CXX_FLAGS_RELEASE}")
endif()

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(PARAMS_FILES )

set(HELPER_FILES
        utils/byteswap.h utils/helper.cpp utils/helper.h utils/testdata.h)

add_executable(bench-dbpit bench/bench-dbpit.cpp ${HELPER_FILES} ${PARAMS_FILES})