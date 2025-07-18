cmake_minimum_required(VERSION 4.0)
project(mpiTest C)

set(CMAKE_C_STANDARD 11)

find_package(MPI REQUIRED)

add_executable(mpiTest
        mpiSample.c)

target_link_libraries(mpiTest PRIVATE MPI::MPI_C)