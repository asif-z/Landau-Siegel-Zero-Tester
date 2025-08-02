//
// Created by kris on 8/2/25.
//

#include "../src/compute.h"
#include "../src/primes.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <flint/long_extras.h>

const slong err = 10;

void test_rhs()
{
    double c = .1;
    double r = 1.45/(10 * log(10));
    printf("%f\n",c);
    printf("%f\n",r);

    arb_t a, b;
    arb_init(a);
    arb_init(b);
    arb_set_d(b,1);
    printf("%d\n",arb_overlaps(a, b));


    compute_config config;
    init_variables(&config);
    for (int q = -1000;  q <= 1000; q += 1)
    {
        if (q == 0 || q== 1) continue;
        double logq = log(abs(q));
        arb_t rhs;
        compute_rhs(&config, q, rhs);
        double rhs_a = (c/r) * (1.0/(r* logq + c)) + (r+ c/logq)/(pow(r+7.0/8.0,2.0))+1.4894 + 0.228774 * logq;
        arb_printd(rhs, 15);
        printf(",%lf,", rhs_a);
        arb_t rhs_ac;
        arb_init(rhs_ac);
        arb_set_d(rhs_ac,rhs_a);
        printf("%d\n",arb_overlaps(rhs, rhs_ac));
        assert(arb_overlaps(rhs, rhs_ac) != 0);
    }
}

void test_compute()
{
    double r = 1.45/(10 * log(10));
    double sigma = 1+r;
    primeiter primes;
    primeiter_init(&primes, "input/primes.txt", 50000);
    long n = 50;
    compute_config config;
    init_variables(&config);
    double sum_a = 0.0;
    for (int q = -100;  q <= 100; q += 1)
    {
        arb_t sum;
        compute_first_n(sum, &config, q, n);
        printf("%d:",q);
        arb_printd(sum, 15);
        sum_a = 0.0;
        set_index(&primes, 0);
        while (primes.index<n)
        {
            double zeta_term = log(primes.cur_prime) *((pow(primes.cur_prime, -sigma))/(1.0-pow(primes.cur_prime, -sigma)));
            double chi = z_kronecker(q, primes.cur_prime);
            double chi_term = log(primes.cur_prime) *((chi * pow(primes.cur_prime, -sigma))/(1.0-chi * pow(primes.cur_prime, -sigma)));
            double term = zeta_term + chi_term;
            sum_a += term;
            get_next_prime(&primes);
        }
        printf(",%lf\n",sum_a);
        arb_t sum_ac;
        arb_init(sum_ac);
        arb_set_d(sum_ac,sum_a);
        arb_add_error_2exp_si(sum_ac, err);
        printf("%d\n",arb_overlaps(sum, sum_ac));
        assert(arb_overlaps(sum, sum_ac) != 0);
    }
}

int main(int argc, char* argv[])
{
    // test_rhs();
    test_compute();
}
