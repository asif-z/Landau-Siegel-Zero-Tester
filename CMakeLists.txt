cmake_minimum_required(VERSION 3.20)
project(mpiTest C)

set(CMAKE_C_STANDARD 11)

find_package(MPI REQUIRED)

enable_testing()

add_subdirectory(test)

add_subdirectory(src)

configure_file(${CMAKE_SOURCE_DIR}/input/primes.txt ${CMAKE_BINARY_DIR}/primes.txt COPYONLY)

configure_file(${CMAKE_SOURCE_DIR}/input/chi.txt ${CMAKE_BINARY_DIR}/chi.txt COPYONLY)

configure_file(${CMAKE_SOURCE_DIR}/input/primes.txt ${CMAKE_BINARY_DIR}/test/primes.txt COPYONLY)

configure_file(${CMAKE_SOURCE_DIR}/input/chi.txt ${CMAKE_BINARY_DIR}/test/chi.txt COPYONLY)

