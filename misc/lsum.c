#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <flint/arb.h>
#include <flint/long_extras.h>
#include "../src/presets.h"
#include "../src/compute.h"
#include <time.h>
#include<unistd.h>
//preset for which value of lambda to use
enum Preset preset = largeX1;

//Computes the zeta and L terms for a small list of d values. Used for debugging and checking how badly some moduli fail

//read a list of primes
int read_primes(long lenPrime, long* primes)
{
    FILE* file;
    int count = 0;
    file = fopen("input/primes.txt", "r");
    if (file == NULL)
    {
        perror("Error opening file");
        return 1;
    }

    // Read each long from the file
    while (fscanf(file, "%ld", &primes[count]) == 1 && count < lenPrime)
    {
        count++;
    }

    fclose(file);
    return 0;
}

int main(int argc, char** argv)
{
    arb_t a;
    arb_init(a);
    arb_log_ui(a,4294967289, 70);
    //Constants

    //len of prime
    const long lenPrime = 5000000;

    //precision set up
    long prec = 70;

    //init var
    arb_t lambda;
    arb_t phi;
    arb_t E;

    arb_init(lambda);
    arb_init(phi);
    arb_init(E);
    initializeLambda(preset, lambda, phi, E, prec);

    long* primes = (long*)malloc(lenPrime * sizeof(long));

    read_primes(lenPrime, primes);

    // setting up variables
    arb_t logd;
    arb_t sigma;
    arb_t sum;
    arb_t logp;
    arb_t con;
    arb_t p;
    arb_t psigma;
    arb_t zeta_term;
    arb_t one; // equal to 1
    arb_init(one);
    arb_set_ui(one, 1);
    arb_t temp1; //temp variables for calculations
    arb_t temp2, temp3, temp4, temp5;
    arb_t l_term;
    arb_t term;
    arb_t c;
    arb_t constant1, constant2, constant3;
    arb_t rhs;
    arb_t r;
    arb_t rhs_term_2;
    arb_t top, bottom;
    arb_t div78; //equal to 7/8
    arb_init(div78);
    arb_set_str(div78, "0.875", prec);

    // sets sigma, r
    arb_init(sigma);
    arb_init(con);
    arb_init(temp1);
    arb_init(r);
    arb_log_ui(con, 10000000000, prec);
    arb_div(r, lambda, con, prec);
    arb_add(sigma, one, r, prec);

    // Set up c
    arb_init(c);
    arb_set_str(c, "0.1", prec);

    long ds[] = {-15,-11,-8,-7,-4,-3 ,5 ,8 ,12 ,13 ,17};
    //loop through all d
    for (int i=0; i<sizeof(ds)/sizeof(long); i++)
    {
        long d = ds[i];
        // Calculating log(d)
        arb_init(logd);
        long absd;
        if (d<0) {
            absd = -d;
        }
        else {
            absd = d;
        }
        arb_log_ui(logd, absd, prec);

        //calculate rhs
        arb_init(rhs);
        arb_init(temp3);
        arb_init(temp4);
        arb_init(temp5);
        arb_init(top);
        arb_init(bottom);
        arb_init(rhs_term_2);

        arb_div(temp3, c, r, prec);
        arb_mul(temp4, r, logd, prec);
        arb_add(temp4, temp4, c, prec);
        arb_div(rhs, temp3, temp4, prec);

        arb_mul(temp5, phi, logd, prec);
        arb_add(rhs, rhs, temp5, prec);
        arb_add(rhs, rhs, E, prec);

        arb_div(top, c, logd, prec);
        arb_add(top, r, top, prec);
        arb_add(bottom, r, div78, prec);
        arb_mul(bottom, bottom, bottom, prec);
        arb_div(rhs_term_2, top, bottom, prec);
        arb_add(rhs, rhs, rhs_term_2, prec);

        //calculate the partial sum
        arb_init(sum);
        arb_t zsum;
        arb_t lsum;
        arb_init(zsum);
        arb_init(lsum);

        for (long primeIndex = 0; primeIndex < lenPrime; primeIndex++)
        {
            long prime = primes[primeIndex];
            int chi = z_kronecker(d, prime); //the value of chi(prime)

            // Calculating log(p)
            arb_init(logp);
            arb_log_ui(logp, prime, prec);

            // Calculating p^sigma
            arb_init(p);
            arb_set_ui(p, prime);
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

            arb_mul(temp1, zeta_term, logp, prec);
            arb_add(zsum, zsum, temp1, prec);
            arb_mul(temp1, l_term, logp, prec);
            arb_add(lsum, lsum, temp1, prec);
        }

        printf("%ld: zsum: ", d);
        printf(arb_get_str(zsum, 5, 0));
        printf(", lsum: ");
        printf(arb_get_str(lsum, 5, 0));
        printf(", rhs: ");
        printf(arb_get_str(rhs, 5, 0));
        printf(", total: ");
        printf(arb_get_str(sum, 5, 0));
        printf("\n");
    }

    return 0;
}
