// primelist.h
#ifndef PRIMES_H
#define PRIMES_H

typedef struct primeiter {
    long* arr;
    long index;
    long cur_prime;
    long size;
} primeiter;

int primeiter_init(primeiter* primes, const char* filename, long lenPrime);

long get_next_prime(primeiter* primes);

long get_prime_at(primeiter* primes, long i);

void set_index(primeiter* primes, long i);

#endif // PRIMES_H