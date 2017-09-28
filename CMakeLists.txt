cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_COMPILER g++-4.8)
set(CMAKE_C_COMPILER g++-4.8)


project(SOP_KMeans CXX)

add_definitions(-DVERSION=\"16.0.671\" -DMAKING_DSO -D_GNU_SOURCE -DLINUX -DAMD64 -m64 -fPIC -DSIZEOF_VOID_P=8 -DFBX_ENABLED=1 -DOPENCL_ENABLED=1 -DOPENVDB_ENABLED=1 -DSESI_LITTLE_ENDIAN -DENABLE_THREADS -DUSE_PTHREADS -D_REENTRANT -D_FILE_OFFSET_BITS=64 -c -DGCC4 -DGCC3 -Wno-deprecated -std=c++11 -isystem /opt/hfs16.0/toolkit/include -Wall -W -Wno-parentheses -Wno-sign-compare -Wno-reorder -Wno-uninitialized -Wunused -Wno-unused-parameter -Wno-unused-local-typedefs -fno-strict-aliasing)
add_definitions(-O3)

LINK_DIRECTORIES(/opt/hfs16.0/dsolib)
INCLUDE_DIRECTORIES(/opt/hfs16.0/toolkit/include)

add_library(${PROJECT_NAME} SHARED "src/SOP_KMean.cpp")
