#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <flint/arb.h>
#include <flint/long_extras.h>
#include <time.h>
#include<unistd.h>

//Computes the zeta and L terms for a small list of q values. Used for debugging and checking how badly some moduli fail


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
    long prec = 70;

    //set up length to calculate
    long qMax = 100000000;

    arb_t lambda;
    arb_t phi;
    arb_t O1; //the O(1)+O(loglog(q)) term
    arb_init(lambda);
    arb_init(phi);
    arb_init(O1);

    //presets:

    //small X (X~2000)
    arb_set_str(lambda, "2.14", prec);
    arb_set_str(phi, "0.219697", prec);
    arb_set_str(O1, "1.33613", prec);

    //large X (X~50000)
    // arb_set_str(lambda, "1.6", prec);
    // arb_set_str(phi, "0.22675", prec);
    // arb_set_str(O1, "1.38794", prec);


    clock_t start, end;
    start = clock(); // start timer on rank 0 only

    long* primes = (long*)malloc(lenPrime * sizeof(long));

    read_primes(lenPrime, primes);

    // setting up variables
    arb_t logq;
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

    long long qs[] = {44701733, -90071083, 86336333, -3, -99989443, -99868267, -98490043, 79188317};
    //loop through all q
    for (int i=0; i<sizeof(qs)/sizeof(long long); i++)
    {
        long long q = qs[i];
        // Calculating log(q)
        arb_init(logq);
        arb_log_ui(logq, abs(q), prec);

        //calculate rhs
        arb_init(rhs);
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

        //calculate the partial sum
        arb_init(sum);
        arb_t zsum;
        arb_t lsum;
        arb_init(zsum);
        arb_init(lsum);

        // loop over all primes < q^alpha
        long prime = primes[0];
        int primeIndex = 0;
        while (primeIndex < 300000)
        {
            int chi = z_kronecker(q, prime); //the value of chi(prime)

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

            prime = primes[primeIndex++];
        }

        printf("%d: zsum: ", q);
        printf(arb_get_str(zsum, 5, 0));
        printf(", lsum: ");
        printf(arb_get_str(lsum, 5, 0));
        printf(", rhs: ");
        printf(arb_get_str(rhs, 5, 0));
        printf(", total: ");
        printf(arb_get_str(sum, 5, 0));
        printf("\n");

        // char *output = (char *) malloc(50 * sizeof(char));
        // output = arb_get_str(rhs, 40, 0);
        // printf("%ld,%s\n",q,output);
    }


    end = clock();
    printf("Elapsed time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);

    return 0;
}
