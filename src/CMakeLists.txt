
add_executable(mpiTest
#        main.c
#        main_test_ind.c
        main_by_file.c
        primes.c
        buffered_chi.c
        compute.c
)

target_include_directories(mpiTest PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(mpiTest PRIVATE flint mpfr gmp m)

target_link_libraries(mpiTest PRIVATE MPI::MPI_C)


