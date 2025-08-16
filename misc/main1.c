#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <flint/arb.h>
#include <flint/ulong_extras.h>
#include <unistd.h>
#include <flint/long_extras.h>
#include "../src/primes.h"
#include "../src/buffered_chi.h"

#define JOB_TAG 1
#define STOP_TAG 2

//total number of primes precomputed
#define lenPrime 500000
//FLINT precision to use
#define prec 50
//max number of primes to add to the sum before truncating
#define primeBd 100000
//maximum modulus used
#define qMax 100
//dimensions of the Kronecker symbol array
#define rows 10000
#define cols 104729

slong const step = 100000;

//Global variables for buffers
buffered_chi chi_value;
primeiter primes;

// setting up gloabl variables
arb_t c; // zero-free region constant
arb_t sigma;
arb_t phi;
arb_t O1; //the O(1)+O(loglog(q)) term
arb_t r;
arb_t div78;
arb_t one;

int init_variables()
{
    if (primeiter_init(&primes, "input/primes.txt",lenPrime) != 0)
    {
        return 1;
    }
    if (chi_init(&chi_value, rows, cols, "input/chi.txt") != 0)
    {
        return 2;
    }

    //Set up zero-free region
    arb_init(c);
    arb_set_str(c, "0.1", prec);

    //init var
    arb_t lambda;
    arb_t logQ;

    arb_init(lambda);
    arb_init(phi);
    arb_init(O1);

    arb_init(one);
    arb_set_ui(one, 1);

    arb_init(div78);
    arb_set_str(div78, "0.875", prec);

    //presets:

    //small X (pi(X)~10000)
    arb_set_str(lambda, "1.4526", prec);
    arb_set_str(phi, "0.22910585", prec);
    arb_set_str(O1, "1.0627878", prec);

    //large X (pi(X)~ 540,000)
    // arb_set_str(lambda, "1.1122", prec);
    // arb_set_str(phi, "0.23344604", prec);
    // arb_set_str(O1, "1.09531399", prec);

    // sets sigma, r
    arb_init(sigma);
    arb_init(logQ);
    arb_init(r);
    arb_log_ui(logQ, 10000000000, prec);
    arb_div(r, lambda, logQ, prec);
    arb_add(sigma, one, r, prec);

    //
    arb_clear(lambda);
    arb_clear(logQ);
    return 0;
}

void compute_rhs(slong q, arb_t rhs)
{
    //init var
    arb_init(rhs);
    arb_t logq;
    arb_t temp3, temp4, temp5, top, bottom, rhs_term_2;
    // Calculating log(q)
    arb_init(logq);
    arb_log_ui(logq, abs(q), prec);

    //calculate rhs
    arb_init(temp3);
    arb_init(temp4);
    arb_init(temp5);
    arb_init(top);
    arb_init(bottom);
    arb_init(rhs_term_2);

    arb_div(temp3, c, r, prec);
    arb_mul(temp4, r, logq, prec);
    arb_add(temp4, temp4, c, prec);
    arb_div(rhs, temp3, temp4, prec);

    arb_mul(temp5, phi, logq, prec);
    arb_add(rhs, rhs, temp5, prec);
    arb_add(rhs, rhs, O1, prec);

    arb_div(top, c, logq, prec);
    arb_add(top, r, top, prec);
    arb_add(bottom, r, div78, prec);
    arb_mul(bottom, bottom, bottom, prec);
    arb_div(rhs_term_2, top, bottom, prec);
    arb_add(rhs, rhs, rhs_term_2, prec);

    //free up space
    arb_clear(logq);
    arb_clear(temp3);
    arb_clear(temp4);
    arb_clear(temp5);
    arb_clear(rhs_term_2);
    arb_clear(bottom);
    arb_clear(top);
}

slong compute(slong q)
{
    arb_t rhs;
    compute_rhs(q, rhs);

    // printf("rhs at %ld:", q);
    // arb_printd(rhs,15);       // prints in standard interval notation
    // printf("\n");

    // loop over primes until we exceed primeBd or the inequality is violated

    //calculate the partial sum
    arb_t sum;
    arb_t logp;
    arb_t p;
    arb_t psigma;
    arb_t zeta_term;
    arb_t temp1;
    arb_t l_term;
    arb_t temp2;
    arb_t term;

    arb_init(sum);

    set_index(&primes, 0);

    while (primes.index < primeBd)
    {
        // compute Kronecker symbol
        int chi = chi_val(&chi_value, q, primes.cur_prime, primes.index);

        // Calculating log(p)
        arb_init(logp);
        arb_log_ui(logp, primes.cur_prime, prec);

        // Calculating p^sigma
        arb_init(p);
        arb_set_ui(p, primes.cur_prime);
        arb_init(psigma);
        arb_pow(psigma, p, sigma, prec);

        // The infinite sum from the p terms for zeta is 1/(p^sigma -1)
        arb_init(zeta_term);
        arb_init(temp1);
        arb_sub(temp1, psigma, one, prec);
        arb_inv(zeta_term, temp1, prec);

        // The infinite sum from the p terms for L is chi(p)/(p^sigma -chi(p))
        arb_init(l_term);
        if (chi == 1)
        {
            arb_set(l_term, zeta_term);
        }
        else if (chi == -1)
        {
            arb_init(temp1);
            arb_init(temp2);
            arb_add(temp1, psigma, one, prec);
            arb_inv(temp2, temp1, prec);
            arb_neg(l_term, temp2);
        }

        // add log(p)*(zeta_term + l_term) to the sum
        arb_init(term);
        arb_init(temp1);
        arb_add(term, zeta_term, l_term, prec);
        arb_mul(temp1, term, logp, prec);
        arb_add(sum, sum, temp1, prec);


        if (primes.index - 50 == 0)
        {
            printf("sum for %d:", q);
            arb_printd(sum, 15); // prints in standard interval notation
            printf("\n");
        }

        get_next_prime(&primes);

        if (primes.index % 50 == 0 && arb_gt(sum, rhs) == 1)
        {
            arb_clear(sum);
            arb_clear(logp);
            arb_clear(temp1);
            arb_clear(temp2);
            arb_clear(term);
            arb_clear(zeta_term);
            arb_clear(l_term);
            arb_clear(p);
            arb_clear(psigma);
            return primes.index;
        }
    }
    arb_clear(sum);
    arb_clear(logp);
    arb_clear(temp1);
    arb_clear(temp2);
    arb_clear(term);
    arb_clear(zeta_term);
    arb_clear(l_term);
    arb_clear(p);
    arb_clear(psigma);
    return -1;
}

int main(int argc, char* argv[])
{
    init_variables();
    compute(-8);
}
