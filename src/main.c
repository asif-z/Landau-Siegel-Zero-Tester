#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>
#include <mpi.h>
#include <stdbool.h>
#include "compute.h"

#define JOB_TAG 1
#define STOP_TAG 2

#define FOLDERNAME_LEN 64

//maximum modulus used
#define qMax 100
slong const step = 100000;

compute_config compute_c;

bool is_valid_q(slong q)
{
    return q != 0 && q != 1 && (q % 4 == 0 || q % 4 == 1 || q % 4 == -3);
}

int master_run(int size)
{
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

    return 0;
}

int worker_run(int rank, char foldername[64])
{
    MPI_Status status;
    slong cur;

    double start = MPI_Wtime();

    printf("Worker %d started\n", rank);
    int code = init_variables(&compute_c);
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

        for (slong q = cur; q < cur + step && q < qMax; q++)
        {
            if (is_valid_q(q))
            {
                slong result = compute(&compute_c, q);
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
