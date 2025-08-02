//
// Created by kris on 8/1/25.
//

#ifndef COMPUTE_H
#define COMPUTE_H

#include <stdlib.h>
#include <flint/arb.h>
#include "primes.h"
#include "buffered_chi.h"

//total number of primes precomputed
#define lenPrime 500000
//FLINT precision to use
#define prec 50
//max number of primes to add to the sum before truncating
#define primeBd 100000
//dimensions of the Kronecker symbol array
#define rows 10000
#define cols 104729

typedef struct compute_config
{
    //Global variables for buffers
    buffered_chi chi_value;
    primeiter primes;

    // setting up gloabl variables
    arb_t c; // zero-free region constant
    arb_t sigma; //1+r
    arb_t phi;
    arb_t O1; //the O(1)+O(loglog(q)) term
    arb_t r;
    arb_t div78; //constant 7/8
    arb_t one; //constant 1
}compute_config;

int init_variables(compute_config *compute_c);

void compute_rhs(compute_config *compute_c,slong q, arb_t rhs);

slong compute(compute_config *compute_c, slong q);

void compute_first_n(arb_t sum, compute_config *compute_c, slong q, slong n);

#endif //COMPUTE_H