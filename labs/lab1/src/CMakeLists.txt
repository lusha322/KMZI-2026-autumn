cmake_minimum_required(VERSION 3.16)
project(aes_gcm_variant1 CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(crypto STATIC
    crypto/gf256.cpp
    crypto/aes128.cpp
    crypto/gcm.cpp
)
target_include_directories(crypto PUBLIC crypto)

add_executable(lab1
    main.cpp
    demo/demo.cpp
    test/test_gf256.cpp
    test/test_aes.cpp
    test/test_gcm.cpp
)
target_link_libraries(lab1 PRIVATE crypto)

enable_testing()
add_test(NAME all_tests COMMAND lab1)