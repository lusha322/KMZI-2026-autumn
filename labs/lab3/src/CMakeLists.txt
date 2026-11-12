cmake_minimum_required(VERSION 3.10)
project(lab3_stb34101_45_variant1 CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# GMP используется только для арифметики больших чисел по модулю q (скаляры подписи)
find_path(GMP_INCLUDE_DIR gmp.h HINTS /opt/homebrew/include /usr/local/include /usr/include)
find_library(GMP_LIBRARY gmp HINTS /opt/homebrew/lib /usr/local/lib /usr/lib)
include_directories(${GMP_INCLUDE_DIR})

add_executable(main
    ./ec_algebra.cpp
    ./ct_secure.cpp
    ./protocols.cpp
    ./main.cpp
)
target_link_libraries(main ${GMP_LIBRARY})
