cmake_minimum_required(VERSION 3.13)
project(ntoolkit VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

set(MY_LIB /mnt/samba/lib)
#include_directories(.
#        ./include/win64
#        /usr/include
#        /usr/local/include)

#link_directories(./lib/win64/tbb-vc14
#        ./bin/win64/tbb-vc14)

set(BOOST_ROOT ${MY_LIB}/boost_1_71_0)

#include_directories(src)
add_subdirectory(src)

add_compile_options(
        -O3
)

add_executable(
        ntoolkit
        lib.cpp
        lib.h
        main.cpp
        model.cpp
        model.h
        neuron.cpp
        neuron.h network.cpp network.h)

# dynamic link?
# link_libraries(tbb.dll)
# static link?
#link_libraries(tbb.lib)
