#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>
#include <mpi.h>
#include <stdbool.h>
#include "compute.h"
#include "presets.h"

/*
 * Uses clusters to verify that L(s,\chi_d) does not have a Landau-Siegel zero for |d|<=qMax
 */

#define JOB_TAG 1
#define STOP_TAG 2

#define FOLDERNAME_LEN 64

//number of bits of precision to use
#define prec0 50
//max number of primes to add to the sum before truncating
#define N00 10000
//number of primes to add to the sum before verifying if it is violated
#define checkDistance0 50
//dimensions of the Kronecker symbol array
#define rows 10000
#define cols 104729
//maximum modulus used
#define qMax 10000
//range of moduli a cluster will compute before requesting more from the master process
#define step 1000
//preset for which value of lambda to use
enum Preset preset = smallX1;

compute_config compute_c0;

int init_variables(compute_config* compute_c)
{
    compute_c->prec = prec0;
    compute_c->N0 = N00;
    compute_c->checkDistance = checkDistance0;

    if (primeiter_init(&(compute_c->primes), "input/primes.txt",N00) != 0)
    {
        return 1;
    }
    if (chi_init(&(compute_c->chi_value), rows, cols, "input/chi.txt") != 0)
    {
        return 2;
    }

    //Set up zero-free region
    arb_init(compute_c->c);
    arb_set_str(compute_c->c, "0.1", prec0);

    //init var
    arb_t lambda;
    arb_t logQ;

    arb_init(lambda);
    arb_init(compute_c->phi);
    arb_init(compute_c->E);
    initializeLambda(preset, lambda, compute_c->phi, compute_c->E, prec0);

    arb_init(compute_c->div78);
    arb_set_str(compute_c->div78, "0.875", prec0);

    // initializes sigma, r
    arb_init(compute_c->sigma);
    arb_init(logQ);
    arb_init(compute_c->r);
    arb_log_ui(logQ, 10000000000, prec0);
    arb_div(compute_c->r, lambda, logQ, prec0);
    arb_add_ui(compute_c->sigma, compute_c->r, 1, prec0);
    arb_clear(lambda);
    arb_clear(logQ);

    compute_zeta_sum(compute_c);
    return 0;
}

// determines whether d is a fundamental discriminant (without the square-free conditions)
bool is_valid_d(long d)
{
    return d != 0 && d != 1 && (d % 16 == 8 || d % 16 == 12 || d % 16 == -8 || d % 16 == -4 || d % 4 == 1 || d % 4 == -
        3);
}

int master_run(int size)
{
    printf("Master started\n");
    fflush(stdout);
    long cur = -qMax; // current modulus

    MPI_Status status;
    int active_workers = size - 1; // Track active workers

    // Assign initial jobs or stop signals
    for (int i = 1; i < size; i++)
    {
        if (cur < qMax)
        {
            MPI_Send(&cur, 1, MPI_LONG_LONG_INT, i, JOB_TAG, MPI_COMM_WORLD);
            cur += step;
        }
        else
        {
            // Send stop immediately if no jobs left
            MPI_Send(NULL, 0, MPI_LONG_LONG_INT, i, STOP_TAG, MPI_COMM_WORLD);
            active_workers--;
        }
    }

    // Manage active workers
    while (active_workers > 0)
    {
        int dummy;
        // Wait for worker ready signal
        MPI_Recv(&dummy, 1, MPI_INT, MPI_ANY_SOURCE, JOB_TAG, MPI_COMM_WORLD, &status);
        int worker_rank = status.MPI_SOURCE;

        if (cur < qMax)
        {
            // Send next job
            MPI_Send(&cur, 1, MPI_LONG_LONG_INT, worker_rank, JOB_TAG, MPI_COMM_WORLD);
            cur += step;
        }
        else
        {
            // Send stop and deactivate worker
            MPI_Send(NULL, 0, MPI_LONG_LONG_INT, worker_rank, STOP_TAG, MPI_COMM_WORLD);
            active_workers--;
        }
    }
    return 0;
}

int worker_run(int rank, char foldername[64])
{
    MPI_Status status;
    long cur;

    double start = MPI_Wtime();

    printf("Worker %d started\n", rank);
    int code = init_variables(&compute_c0);
    if (code != 0)
    {
        MPI_Abort(MPI_COMM_WORLD, code);
        return 1;
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

        for (long d = cur; d < cur + step && d < qMax; d++)
        {
            if (is_valid_d(d))
            {
                long result = compute(&compute_c0, d);
                if (result < 0)
                {
                    fprintf(outfile, "%ld,fail,%ld\n", d, result);
                }
                else
                {
                    fprintf(outfile, "%ld,pass,%ld\n", d, result);
                }
            }
        }
        fflush(outfile);

        // Signal master that this worker is ready
        int ready = 1;
        MPI_Send(&ready, 1, MPI_INT, 0, JOB_TAG, MPI_COMM_WORLD);
    }

    fclose(outfile);
    return 0;
}

void create_and_broadcast_folder(char foldername[FOLDERNAME_LEN], int rank)
{
    if (rank == 0)
    {
        // Create timestamped folder name
        time_t t = time(NULL);
        struct tm tm = *localtime(&t);
        snprintf(foldername, FOLDERNAME_LEN, "run_%04d%02d%02d_%02d%02d%02d",
                 tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday,
                 tm.tm_hour, tm.tm_min, tm.tm_sec);

        // Create the folder
        if (mkdir(foldername, 0755) != 0)
        {
            perror("mkdir failed on rank 0");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // Broadcast the folder name to all other ranks
    MPI_Bcast(foldername, FOLDERNAME_LEN, MPI_CHAR, 0, MPI_COMM_WORLD);
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

    char foldername[FOLDERNAME_LEN];

    // Rank 0 creates folder, others receive name
    create_and_broadcast_folder(foldername, rank);

    MPI_Barrier(MPI_COMM_WORLD);
    double start;
    if (rank == 0)
    {
        start = MPI_Wtime();
    }

    if (rank == 0)
    {
        master_run(size);
    }
    else
    {
        // Worker process
        worker_run(rank, foldername);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0)
    {
        double total_time = MPI_Wtime() - start;
        printf("Total Time: %.3f seconds\n", total_time);
    }

    MPI_Finalize();
    return 0;
}
