//
// Created by kris on 8/1/25.
//
#include <stdio.h>
#include <stdlib.h>
#include "primes.h"

//reads a precomputed list of primes
int primeiter_init(primeiter* primes, const char* filename, long size)
{
    primes->arr = (long*)malloc(size * sizeof(long));

    if (primes->arr == NULL)
    {
        perror("Memory allocation failed");
        return 1;
    }

    FILE* file;
    int count = 0;
    file = fopen(filename, "r");
    if (file == NULL)
    {
        perror("Error opening file");
        return 1;
    }

    // Read each long from the file
    while (count < size && fscanf(file, "%ld", &primes->arr[count]) == 1)
    {
        count++;
    }

    fclose(file);

    primes->index = 0;
    primes->cur_prime = 2;
    primes->size = size;
    return 0;
}

long get_next_prime(primeiter* primes)
{
    if (primes->index == primes->size - 1)
    {
        return -1;
    }
    primes->index++;
    primes->cur_prime = primes->arr[primes->index];
    return primes->index;
}

long get_prime_at(primeiter* primes, long i)
{
    return primes->arr[i];
}

void set_index(primeiter* primes, long i)
{
    primes->index = i;
    primes->cur_prime = primes->arr[i];
}
