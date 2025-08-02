
#include <stdio.h>
#include <flint/ulong_extras.h>
#include "../src/primes.h"
#include <assert.h>

int main(int argc, char* argv[])
{
    primeiter primes;
    primeiter_init(&primes, "input/primes.txt", 50000);
    n_primes_t iter;
    n_primes_init(iter);
    do
    {
        assert(n_primes_next(iter) == primes.cur_prime);
        printf("%ld\n",primes.cur_prime);
    }
    while (get_next_prime(&primes) != -1);
}
