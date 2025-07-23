#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <flint/arb.h>
#include <flint/long_extras.h>
#include <time.h>
#include<unistd.h>

//read a list of primes
int read_primes(long lenPrime, long* primes)
{
    FILE* file;
    int count = 0;
    file = fopen("primes.txt", "r");
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
    //Constants

    //len of prime
    const long lenPrime = 500000;

    //precision set up
    long prec = 50;

    //set up lambda
    arb_t lambda;
    arb_init(lambda);
    arb_set_str(lambda, "0.69165", prec);

    long* primes = (long*)malloc(lenPrime * sizeof(long));
    read_primes(lenPrime, primes);

    // setting up variables
    arb_t sigma;
    arb_t logp;
    arb_t logQ;
    arb_t p;
    arb_t psigma;
    arb_t zeta_term;
    arb_t one; // equal to 1
    arb_init(one);
    arb_set_ui(one, 1);
    arb_t temp1; //temp variables for calculations
    arb_t temp2, temp3,  temp4, temp5;
    arb_t r;
    arb_t zsum;

    // sets sigma, r
    arb_init(sigma);
    arb_init(logQ);
    arb_init(temp1);
    arb_init(r);
    arb_log_ui(logQ, 10000000000, prec);
    arb_div(r, lambda, logQ, prec);
    arb_add(sigma, one, r, prec);

            
    arb_init(zsum);
    long prime = primes[0];
    int primeIndex = 0;
    while (primeIndex < 10000)
    {
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

        // add log(p)*(zeta_term + l_term) to the sum
        arb_mul(temp1, zeta_term, logp, prec);
        arb_add(zsum, zsum, temp1, prec);

        prime = primes[primeIndex++];

        if(primeIndex==100 || primeIndex%500==0){
            printf("%d,%d: zsum: ", primeIndex, prime);
            printf(arb_get_str(zsum, 5, 0));
            printf("\n");
        }
    }

    return 0;
}
