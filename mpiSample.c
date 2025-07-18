#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <flint/arb.h>
#include <flint/long_extras.h>
#include <time.h>

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
    const long lenPrime = 50000;

    //precision set up
    long prec = 100;

    //set up alpha
    double alpha = 0.3;

    //set up length to calculate
    long qMax = 1000000;

    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //set up lambda
    arb_t lambda;
    arb_init(lambda);
    arb_set_d(lambda, 1);


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

    printf("Elapsed time after file reading: %.3f seconds\n", (double)(clock() - start) / CLOCKS_PER_SEC);

    // setting up variables
    arb_t logq;
    arb_t sigma;
    arb_t sum;
    arb_t logp;
    arb_t p;
    arb_t psigma;
    arb_t zeta_term;
    arb_t one; // equal to 1
    arb_init(one);
    arb_set_ui(one, 1);
    arb_t temp1; //temp variables for calculations
    arb_t temp2;
    arb_t l_term;
    arb_t term;

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
        if ((q + qMax) % size != rank) continue; // skip q not assigned to this rank

        if (q % 100000 == 0)
        {
            printf("Elapsed time when q=%d: %.3f seconds\n", q, (double)(clock() - start) / CLOCKS_PER_SEC);
        }

        if (q % 4 == 0 || q % 4 == 1)
        {
            //number of terms to compute
            long len = pow(abs(q), alpha);

            // sets sigma= 1 + lambda/log(q)
            arb_init(sigma);
            arb_init(logq);
            arb_init(temp1);
            arb_log_ui(logq, abs(q), prec);
            arb_div(temp1, lambda, logq, prec);
            arb_add(sigma, one, temp1, prec);

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
                else
                {
                }

                // add log(p)*(zeta_term + l_term) to the sum
                arb_init(term);
                arb_init(temp1);
                arb_add(term, zeta_term, l_term, prec);
                arb_mul(temp1, term, logp, prec);
                arb_add(sum, sum, temp1, prec);

                prime = primes[primeIndex++];
            }

            char* output = (char*)malloc(50 * sizeof(char));
            output = arb_get_str(sum, 40, 0);
            fprintf(outfile, "%ld,%s\n", q, output);
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
