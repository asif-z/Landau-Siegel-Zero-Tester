#include <stdio.h>

#include "../src/buffered_chi.h"
#include "../src/primes.h"
#include <assert.h>
#include <flint/long_extras.h>


#define rows 10000
#define cols 104729

int main(int argc, char* argv[])
{
    buffered_chi chi;
    chi_init(&chi, rows, cols, "input/chi.txt");

    primeiter primes;
    primeiter_init(&primes, "input/primes.txt", 50000);
    for (int q = -100; q < 100; q++)
    {
        do
        {
            assert(chi_val(&chi, q, primes.cur_prime, primes.index) == z_kronecker(q, primes.cur_prime));
        }
        while (get_next_prime(&primes) != -1);
    }
}
