// primelist.h
#ifndef PRIMES_H
#define PRIMES_H

typedef struct primeiter {
    long* arr; //array of primes
    long index; //current index
    long cur_prime; //current prime
    long size; //size of the array
} primeiter;

int primeiter_init(primeiter* primes, const char* filename, long size);

long get_next_prime(primeiter* primes);

long get_prime_at(primeiter* primes, long i);

void set_index(primeiter* primes, long i);

#endif // PRIMES_H