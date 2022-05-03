cmake_minimum_required(VERSION 3.20)
project(TBlib)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "-mavx -mavx2 -mfma -msse -msse2 -msse3")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${C_FLAGS} ${CXX_FLAGS}")

include_directories(include include/tbl include/tbl/prng include/tbl/helper WinNTL-11_4_3/include include/tbl/layer3)

file(GLOB TBLIB_SRCS lib/*.cpp)
file(GLOB NTL_SRCS WinNTL-11_4_3/src/*.cpp)

add_library(TBlib_extra SHARED ${TBLIB_SRCS})
add_library(NTL STATIC ${NTL_SRCS})
add_executable(TBlib main.cpp utils/bench.h utils/bench.cpp)
target_link_libraries(TBlib NTL TBlib_extra)

if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(C_FLAGS "${C_FLAGS} -O2 -fvisibility=hidden -fvisibility-inlines-hidden")
    set(LINKER_FLAGS "${LINKER_FLAGS} -Wl,-exclude-libs,ALL -Wl,--gc-sections -Wl,--strip-all")
else ()
    set(C_FLAGS "${C_FLAGS} -O0")
endif ()