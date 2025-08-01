
add_executable(mpiTest
        main.c
        primes.c
        buffered_chi.c
)

target_include_directories(mpiTest PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(mpiTest PRIVATE flint mpfr gmp m)

target_link_libraries(mpiTest PRIVATE MPI::MPI_C)


