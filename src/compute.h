//
// Created by kris on 8/1/25.
//

#ifndef COMPUTE_H
#define COMPUTE_H

#include <stdlib.h>
#include <flint/arb.h>
#include "primes.h"
#include "buffered_chi.h"

typedef struct compute_config
{
    long N0; //the maximum number of primes we will sum up to
    long prec; // number of bits of precision to use in calculations
    long checkDistance; // number of primes to add to the sum before verifying if it is violated

    //Global variables for buffers
    buffered_chi chi_value;
    primeiter primes;
    arb_t* zetaSums;

    // setting up global variables
    arb_t c; // the zero-free region constant
    arb_t sigma; //=1+r
    arb_t phi; //the coefficient on log(q)
    arb_t E; //the O(1)+O(log(1+1/r)) term
    arb_t r;
    arb_t div78; //constant 7/8
}compute_config;

void compute_rhs(compute_config *compute_c, long d, arb_t rhs);

void compute_zeta_sum(compute_config *compute_c);

long compute(compute_config *compute_c, long d);

void compute_first_n(arb_t sum, compute_config *compute_c, long d, long n);

#endif //COMPUTE_H