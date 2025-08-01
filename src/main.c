#include <stdio.h>
#include <stdlib.h>
#include <flint/arb.h>
#include <sys/stat.h>
#include <time.h>
#include <mpi.h>
#include <stdbool.h>
#include "primes.h"
#include "buffered_chi.h"

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
arb_t sigma; //1+r
arb_t phi;
arb_t O1; //the O(1)+O(loglog(q)) term
arb_t r;
arb_t div78; //constant 7/8
arb_t one; //constant 1

int init_variables()
{
    if (primeiter_init(&primes, "input/primes.txt",primeBd) != 0)
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

    //small X (pi(X)~7000)
    arb_set_str(lambda, "1.45", prec);
    arb_set_str(phi, "0.228774", prec);
    arb_set_str(O1, "1.4894", prec);

    //large X (pi(X)~200000)
    // arb_set_str(lambda, "1.3", prec);
    // arb_set_str(phi, "0.23083", prec);
    // arb_set_str(O1, "1.50458", prec);

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

bool is_valid_q(slong q) {
    return q != 0 && q != 1 && (q % 4 == 0 || q % 4 == 1 || q % 4 == -3);
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


        // if (primes.index - 50 == 0)
        // {
        //     printf("sum for %ld:", q);
        //     arb_printd(sum, 15); // prints in standard interval notation
        //     printf("\n");
        // }

        if (get_next_prime(&primes) == -1)
        {
            break;
        }

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

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size < 2)
    {
        printf("2 or more mpi processes are needed");
        MPI_Finalize();
        return 1;
    }

    char foldername[64];
    // Create a unique folder name using timestamp
    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    snprintf(foldername, sizeof(foldername), "run_%04d%02d%02d_%02d%02d%02d",
             tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday,
             tm.tm_hour, tm.tm_min, tm.tm_sec);

    // Create the folder
    mkdir(foldername, 0755);


    if (rank == 0)
    {
        double start = MPI_Wtime();
        printf("Master started\n");
        fflush(stdout);
        slong cur = -qMax;
        MPI_Status status;

        // Assign initial jobs
        for (int i = 1; i < size && cur < qMax; i++)
        {
            MPI_Send(&cur, 1, MPI_LONG_LONG_INT, i, JOB_TAG, MPI_COMM_WORLD);
            cur += step;
        }

        // Continue assigning jobs as workers become available
        while (cur < qMax)
        {
            // Receive a ready signal (just a dummy receive)
            int dummy;
            MPI_Recv(&dummy, 1, MPI_INT, MPI_ANY_SOURCE, JOB_TAG, MPI_COMM_WORLD, &status);
            MPI_Send(&cur, 1, MPI_LONG_LONG_INT, status.MPI_SOURCE, JOB_TAG, MPI_COMM_WORLD);
            cur += step;
        }

        // Send stop signal to all workers
        for (int i = 1; i < size; i++)
        {
            MPI_Send(NULL, 0, MPI_LONG_LONG_INT, i, STOP_TAG, MPI_COMM_WORLD);
        }

        printf("Total Time: %f", MPI_Wtime() - start);
    }
    else
    {
        // Worker process
        MPI_Status status;
        slong cur;

        double start = MPI_Wtime();

        printf("Worker %d started\n", rank);
        int code = init_variables();
        if (code != 0)
        {
            MPI_Abort(MPI_COMM_WORLD, code);
            return code;
        }

        printf("Files loaded at rank %d: %f\n", rank, MPI_Wtime() - start);

        FILE* outfile;
        char filename[256];
        snprintf(filename, sizeof(filename), "%s/output_rank_%d.csv", foldername, rank);
        outfile = fopen(filename, "w");

        if (!outfile)
        {
            fprintf(stderr, "Worker %d: Failed to open output file.\n", rank);
            MPI_Finalize();
            return 1;
        }

        while (1)
        {
            MPI_Recv(&cur, 1, MPI_LONG_LONG_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            if (status.MPI_TAG == STOP_TAG) break;

            // Process job
            printf("%ld at worker %d\n", cur, rank);

            for (slong q = cur; q < cur + step && q < qMax; q++)
            {
                if (is_valid_q(q))
                {
                    slong result = compute(q);
                    if (result < 0)
                    {
                        fprintf(outfile, "%ld,fail,%ld\n", q, result);
                    }
                    else
                    {
                        fprintf(outfile, "%ld,pass,%ld\n", q, result);
                    }
                }
            }
            fflush(outfile);

            // Signal master that this worker is ready
            int ready = 1;
            MPI_Send(&ready, 1, MPI_INT, 0, JOB_TAG, MPI_COMM_WORLD);
        }

        fclose(outfile);
    }

    MPI_Finalize();
    return 0;
}
