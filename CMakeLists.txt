cmake_minimum_required(VERSION 3.10)
project(mpiTest C)

set(CMAKE_C_STANDARD 11)

find_package(MPI REQUIRED)

add_executable(mpiTest
        main.c)

target_link_libraries(mpiTest PRIVATE MPI::MPI_C)

target_link_libraries(mpiTest PRIVATE flint mpfr gmp m)

configure_file(${CMAKE_SOURCE_DIR}/primes.txt ${CMAKE_BINARY_DIR}/primes.txt COPYONLY)