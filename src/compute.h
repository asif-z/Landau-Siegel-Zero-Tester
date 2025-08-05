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
    long primeBd;
    long prec;

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

void compute_rhs(compute_config *compute_c, long q, arb_t rhs);

long compute(compute_config *compute_c, long q);

void compute_first_n(arb_t sum, compute_config *compute_c, long q, long n);

#endif //COMPUTE_H