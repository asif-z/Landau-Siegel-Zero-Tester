#include <mpi.h>
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

    //set up alpha
    double alpha = .6;

    //set up length to calculate
    long qMax = 100000000;

    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //set up lambda
    arb_t lambda;
    arb_init(lambda);
    arb_set_str(lambda, "1.65", prec);


    clock_t start, end;
    if (rank == 0)
    {
        start = clock(); // start timer on rank 0 only
    }

    long* primes = (long*)malloc(lenPrime * sizeof(long));

    if (read_primes(lenPrime, primes) == 1)
    {
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }

    if (rank == 0)
    {
        printf("Elapsed time after file reading: %.3f seconds\n", (double)(clock() - start) / CLOCKS_PER_SEC);
    }

    // setting up variables
    arb_t logq;
    arb_t sigma;
    arb_t sum;
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
    arb_t l_term;
    arb_t term;
    arb_t c;
    arb_t phi;
    arb_t O1;
    arb_t rhs;
    arb_t r;

    // sets sigma, r
    arb_init(sigma);
    arb_init(logQ);
    arb_init(temp1);
    arb_init(r);
    arb_log_ui(logQ, 10000000000, prec);
    arb_div(r, lambda, logQ, prec);
    arb_add(sigma, one, r, prec);

    // Set up c
    arb_init(c);
    arb_set_str(c, "0.01", prec);

    arb_init(phi);
    arb_init(O1);
    arb_set_str(phi, "0.2432", prec);
    arb_set_str(O1, "2.8943", prec);


    // Open separate output file per rank to avoid clashes
    char filename[100];
    snprintf(filename, sizeof(filename), "output_rank_%d.csv", rank);
    FILE* outfile = fopen(filename, "w");
    if (outfile == NULL)
    {
        perror("Error opening file");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    //loop through all q
    for (long q = -qMax; q <= qMax; q++)
    {
        if (q % 100000 == 0 & rank == 0)
        {
            printf("Elapsed time when q=%d: %.3f seconds\n", q, (double)(clock() - start) / CLOCKS_PER_SEC);
        }

        if ((q + qMax) % size != rank) continue; // skip q not assigned to this rank

        if (q==0) continue;

        // Calculating log(q)
        arb_init(logq);
        arb_log_ui(logq, abs(q), prec);

        if (q % 4 == 0 || q % 4 == 1)
        {
            //number of terms to compute
            long len = pow(abs(q), alpha);

            //calculate the partial sum
            arb_init(sum);

            // loop over all primes < q^alpha
            long prime = primes[0];
            int primeIndex = 0;
            while (prime < len)
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

                prime = primes[primeIndex++];
            }

            arb_init(rhs);
            arb_init(temp3);
            arb_init(temp4);
            arb_init(temp5);

            arb_div(temp3, c, r, prec);
            arb_mul(temp4, r, logq, prec);
            arb_add(temp4, temp4, c, prec);
            arb_div(rhs, temp3, temp4, prec);
            arb_mul(temp5, phi, logq, prec);
            arb_add(rhs, rhs, temp5, prec);
            arb_add(rhs, rhs, O1, prec);

            // char *output = (char *) malloc(50 * sizeof(char));
            // output = arb_get_str(rhs, 40, 0);
            // printf("%ld,%s\n",q,output);

            if (arb_gt(sum, rhs) == 1)
            {
                fprintf(outfile, "%ld,pass\n", q);
            }else
            {
                fprintf(outfile, "%ld,fail\n", q);
            }

        }
    }

    fclose(outfile);

    if (rank == 0)
    {
        end = clock();
        printf("Elapsed time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    }

    MPI_Finalize();

    return 0;
}
