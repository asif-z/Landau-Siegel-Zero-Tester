
#include "compute.h"



int init_variables(compute_config *compute_c)
{
    if (primeiter_init(&(compute_c->primes), "input/primes.txt",primeBd) != 0)
    {
        return 1;
    }
    if (chi_init(&(compute_c->chi_value), rows, cols, "input/chi.txt") != 0)
    {
        return 2;
    }

    //Set up zero-free region
    arb_init(compute_c->c);
    arb_set_str(compute_c->c, "0.1", prec);

    //init var
    arb_t lambda;
    arb_t logQ;

    arb_init(lambda);
    arb_init(compute_c->phi);
    arb_init(compute_c->O1);

    arb_init(compute_c->one);
    arb_set_ui(compute_c->one, 1);

    arb_init(compute_c->div78);
    arb_set_str(compute_c->div78, "0.875", prec);

    //presets:

    //small X (pi(X)~7000)
    arb_set_str(lambda, "1.45", prec);
    arb_set_str(compute_c->phi, "0.228774", prec);
    arb_set_str(compute_c->O1, "1.4894", prec);

    //large X (pi(X)~200000)
    // arb_set_str(lambda, "1.3", prec);
    // arb_set_str(phi, "0.23083", prec);
    // arb_set_str(O1, "1.50458", prec);

    // sets sigma, r
    arb_init(compute_c->sigma);
    arb_init(logQ);
    arb_init(compute_c->r);
    arb_log_ui(logQ, 10000000000, prec);
    arb_div(compute_c->r, lambda, logQ, prec);
    arb_add(compute_c->sigma, compute_c->one, compute_c->r, prec);

    //
    arb_clear(lambda);
    arb_clear(logQ);
    return 0;
}

void compute_rhs(compute_config *compute_c,slong q, arb_t rhs)
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

    arb_div(temp3, compute_c->c, compute_c->r, prec);
    arb_mul(temp4, compute_c->r, logq, prec);
    arb_add(temp4, temp4, compute_c->c, prec);
    arb_div(rhs, temp3, temp4, prec);

    arb_mul(temp5, compute_c->phi, logq, prec);
    arb_add(rhs, rhs, temp5, prec);
    arb_add(rhs, rhs, compute_c->O1, prec);

    arb_div(top, compute_c->c, logq, prec);
    arb_add(top, compute_c->r, top, prec);
    arb_add(bottom, compute_c->r, compute_c->div78, prec);
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

slong compute(compute_config *compute_c, slong q)
{
    arb_t rhs;
    compute_rhs(compute_c, q, rhs);

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

    set_index(&(compute_c->primes), 0);

    while (compute_c->primes.index < primeBd)
    {
        // compute Kronecker symbol
        int chi = chi_val(&compute_c->chi_value, q, compute_c->primes.cur_prime, compute_c->primes.index);

        // Calculating log(p)
        arb_init(logp);
        arb_log_ui(logp, compute_c->primes.cur_prime, prec);

        // Calculating p^sigma
        arb_init(p);
        arb_set_ui(p, compute_c->primes.cur_prime);
        arb_init(psigma);
        arb_pow(psigma, p, compute_c->sigma, prec);

        // The infinite sum from the p terms for zeta is 1/(p^sigma -1)
        arb_init(zeta_term);
        arb_init(temp1);
        arb_sub(temp1, psigma, compute_c->one, prec);
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
            arb_add(temp1, psigma, compute_c->one, prec);
            arb_inv(temp2, temp1, prec);
            arb_neg(l_term, temp2);
        }

        // add log(p)*(zeta_term + l_term) to the sum
        arb_init(term);
        arb_init(temp1);
        arb_add(term, zeta_term, l_term, prec);
        arb_mul(temp1, term, logp, prec);
        arb_add(sum, sum, temp1, prec);


        // if (primes.index - 50 == 0)
        // {
        //     printf("sum for %ld:", q);
        //     arb_printd(sum, 15); // prints in standard interval notation
        //     printf("\n");
        // }

        if (get_next_prime(&compute_c->primes) == -1)
        {
            break;
        }

        if (compute_c->primes.index % 50 == 0 && arb_gt(sum, rhs) == 1)
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
            return compute_c->primes.index;
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