cmake_minimum_required(VERSION 3.23)
project(HW2)

set(CMAKE_CXX_STANDARD 14)

include_directories(/usr/local/include)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xclang -fopenmp")
find_library(OPENMP_LIB libomp.dylib /usr/local/lib)

add_executable(HW2
        main.cpp Grid2d.cpp Grid2d.h math_tools.cpp math_tools.h SL_method.cpp SL_method.h hw2.q1.cpp hw2.q2.cpp)

target_link_libraries(HW2 LINK_PUBLIC ${OPENMP_LIB})
