
#include "compute.h"

// calculates the right side of equation (5) in the paper
void compute_rhs(compute_config *compute_c, long d, arb_t rhs)
{
    //init var
    arb_init(rhs);
    arb_t logq;
    arb_t temp3, temp4, temp5, top, bottom, rhs_term_2;
    // Calculating log(|d|)
    arb_init(logq);
    long absd;
    if (d<0) {
        absd = -d;
    }
    else {
        absd = d;
    }
    arb_log_ui(logq, absd, compute_c->prec);

    //calculate rhs
    arb_init(temp3);
    arb_init(temp4);
    arb_init(temp5);
    arb_init(top);
    arb_init(bottom);
    arb_init(rhs_term_2);

    arb_div(temp3, compute_c->c, compute_c->r, compute_c->prec);
    arb_mul(temp4, compute_c->r, logq, compute_c->prec);
    arb_add(temp4, temp4, compute_c->c, compute_c->prec);
    arb_div(rhs, temp3, temp4, compute_c->prec);

    arb_mul(temp5, compute_c->phi, logq, compute_c->prec);
    arb_add(rhs, rhs, temp5, compute_c->prec);
    arb_add(rhs, rhs, compute_c->E, compute_c->prec);

    arb_div(top, compute_c->c, logq, compute_c->prec);
    arb_add(top, compute_c->r, top, compute_c->prec);
    arb_add(bottom, compute_c->r, compute_c->div78, compute_c->prec);
    arb_mul(bottom, bottom, bottom, compute_c->prec);
    arb_div(rhs_term_2, top, bottom, compute_c->prec);
    arb_add(rhs, rhs, rhs_term_2, compute_c->prec);

    //free up space
    arb_clear(logq);
    arb_clear(temp3);
    arb_clear(temp4);
    arb_clear(temp5);
    arb_clear(rhs_term_2);
    arb_clear(bottom);
    arb_clear(top);
}

void compute_zeta_sum(compute_config *compute_c) {
    compute_c->zetaSums = malloc(compute_c->N0/compute_c->checkDistance*sizeof(arb_t));

    arb_t sum;
    arb_t logp;
    arb_t p;
    arb_t psigma;
    arb_t zeta_term;
    arb_t temp1;
    arb_t temp2;
    arb_init(sum);
    arb_init(temp1);
    arb_init(temp2);

    set_index(&(compute_c->primes), 0);

    while (compute_c->primes.index < compute_c->N0)
    {
        // Calculating log(p)
        arb_init(logp);
        arb_log_ui(logp, compute_c->primes.cur_prime, compute_c->prec);

        // Calculating p^sigma
        arb_init(p);
        arb_set_ui(p, compute_c->primes.cur_prime);
        arb_init(psigma);
        arb_pow(psigma, p, compute_c->sigma, compute_c->prec);

        // The infinite sum from the p terms for zeta is 1/(p^sigma -1)
        arb_init(zeta_term);
        arb_init(temp1);
        arb_sub_ui(temp1, psigma, 1, compute_c->prec);
        arb_inv(zeta_term, temp1, compute_c->prec);

        // add log(p)*zeta_term to the sum
        arb_init(temp1);
        arb_mul(temp1, zeta_term, logp, compute_c->prec);
        arb_add(sum, sum, temp1, compute_c->prec);

        // add a value to the array every checkDistance primes
        if (compute_c->primes.index % compute_c->checkDistance == 0)
        {
            arb_init(compute_c->zetaSums[compute_c->primes.index/compute_c->checkDistance]);
            arb_set(compute_c->zetaSums[compute_c->primes.index/compute_c->checkDistance], sum);
        }

        if (get_next_prime(&compute_c->primes) == -1)
        {
            break;
        }
    }
    arb_clear(sum);
    arb_clear(logp);
    arb_clear(temp1);
    arb_clear(temp2);
    arb_clear(zeta_term);
    arb_clear(p);
    arb_clear(psigma);
}

// Calculates the sum in equation (8) until the inequality is violated or the cut-off value is reached
long compute(compute_config *compute_c, long d)
{
    arb_t rhs;
    compute_rhs(compute_c, d, rhs);

    // loop over primes until we exceed N0 or the inequality is violated

    //calculate the partial sum
    arb_t sum;
    arb_t logp;
    arb_t p;
    arb_t psigma;
    arb_t temp1;
    arb_t l_term;
    arb_t temp2;

    arb_init(sum);
    arb_init(temp1);
    arb_init(temp2);

    set_index(&(compute_c->primes), 0);

    while (compute_c->primes.index < compute_c->N0)
    {
        // compute Kronecker symbol
        int chi = chi_val(&compute_c->chi_value, d, compute_c->primes.cur_prime, compute_c->primes.index);

        // Calculating log(p)
        arb_init(logp);
        arb_log_ui(logp, compute_c->primes.cur_prime, compute_c->prec);

        // Calculating p^sigma
        arb_init(p);
        arb_set_ui(p, compute_c->primes.cur_prime);
        arb_init(psigma);
        arb_pow(psigma, p, compute_c->sigma, compute_c->prec);

        // The infinite sum from the p terms for L is chi(p)/(p^sigma -chi(p))
        arb_init(l_term);
        if (chi == 1)
        {
            arb_init(temp1);
            arb_sub_ui(temp1, psigma, 1, compute_c->prec);
            arb_inv(l_term, temp1, compute_c->prec);
        }
        else if (chi == -1)
        {
            arb_init(temp1);
            arb_init(temp2);
            arb_add_ui(temp1, psigma, 1, compute_c->prec);
            arb_inv(temp2, temp1, compute_c->prec);
            arb_neg(l_term, temp2);
        }

        // add log(p)*l_term to the sum
        arb_init(temp1);
        arb_mul(temp1, l_term, logp, compute_c->prec);
        arb_add(sum, sum, temp1, compute_c->prec);

        // check if the inequality is violated every checkDistance primes
        if (compute_c->primes.index % compute_c->checkDistance == 0)
        {
            arb_add(temp1, sum, compute_c->zetaSums[compute_c->primes.index/compute_c->checkDistance], compute_c->prec);
            if (arb_gt(temp1, rhs) == 1) {
                arb_clear(sum);
                arb_clear(logp);
                arb_clear(temp1);
                arb_clear(temp2);
                arb_clear(l_term);
                arb_clear(p);
                arb_clear(psigma);
                return compute_c->primes.index;
            }
        }

        if (get_next_prime(&compute_c->primes) == -1)
        {
            break;
        }
    }
    arb_clear(sum);
    arb_clear(logp);
    arb_clear(temp1);
    arb_clear(temp2);
    arb_clear(l_term);
    arb_clear(p);
    arb_clear(psigma);
    return -1;
}

// used for testing purposes. Calculates the sum in equation (8) for n primes, and then returns the sum
void compute_first_n(arb_t sum, compute_config *compute_c, long d, long n)
{
    arb_t rhs;
    compute_rhs(compute_c, d, rhs);

    // loop over primes until we exceed N0 or the inequality is violated

    //calculate the partial sum
    arb_t logp;
    arb_t p;
    arb_t psigma;
    arb_t zeta_term;
    arb_t temp1;
    arb_t l_term;
    arb_t temp2;
    arb_t term;

    arb_init(sum);
    arb_init(temp1);
    arb_init(temp2);

    set_index(&(compute_c->primes), 0);

    while (compute_c->primes.index < n)
    {
        // compute Kronecker symbol
        int chi = chi_val(&compute_c->chi_value, d, compute_c->primes.cur_prime, compute_c->primes.index);

        // Calculating log(p)
        arb_init(logp);
        arb_log_ui(logp, compute_c->primes.cur_prime, compute_c->prec);

        // Calculating p^sigma
        arb_init(p);
        arb_set_ui(p, compute_c->primes.cur_prime);
        arb_init(psigma);
        arb_pow(psigma, p, compute_c->sigma, compute_c->prec);

        // The infinite sum from the p terms for zeta is 1/(p^sigma -1)
        arb_init(zeta_term);
        arb_init(temp1);
        arb_sub_ui(temp1, psigma, 1, compute_c->prec);
        arb_inv(zeta_term, temp1, compute_c->prec);

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
            arb_add_ui(temp1, psigma, 1, compute_c->prec);
            arb_inv(temp2, temp1, compute_c->prec);
            arb_neg(l_term, temp2);
        }

        // add log(p)*(zeta_term + l_term) to the sum
        arb_init(term);
        arb_init(temp1);
        arb_add(term, zeta_term, l_term, compute_c->prec);
        arb_mul(temp1, term, logp, compute_c->prec);
        arb_add(sum, sum, temp1, compute_c->prec);

        if (get_next_prime(&compute_c->primes) == -1)
        {
            printf("out of bound!");
            break;
        }
    }
    arb_clear(logp);
    arb_clear(temp1);
    arb_clear(temp2);
    arb_clear(term);
    arb_clear(zeta_term);
    arb_clear(l_term);
    arb_clear(p);
    arb_clear(psigma);
}