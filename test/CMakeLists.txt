# test/CMakeLists.txt

add_executable(unitTest
#        primes_test.c
        ../src/primes.c
#        chi_test.c
        ../src/buffered_chi.c
        ../src/compute.c
        ../src/presets.c
        computation_test.c
)

target_include_directories(unitTest PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(unitTest PRIVATE flint mpfr gmp m)

add_test(NAME UnitTest COMMAND unitTest)
