#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "compute.h"
#include "presets.h"

#define MAX_NUMBERS 100000  // Change this as needed

//number of bits of precision to use
#define prec0 50
//max number of primes to add to the sum before truncating
#define N00 5000000
//number of primes to add to the sum before verifying if it is violated
#define checkDistance0 5000
//dimensions of the Kronecker symbol array
#define rows 10000
#define cols 104729
//preset for which value of lambda to use
enum Preset preset = smallX1;

int init_variables(compute_config *compute_c)
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

    // sets sigma, r
    arb_init(compute_c->sigma);
    arb_init(logQ);
    arb_init(compute_c->r);
    arb_log_ui(logQ, 10000000000, prec0);
    arb_div(compute_c->r, lambda, logQ, prec0);
    arb_add_ui(compute_c->sigma, compute_c->r, 1, prec0);

    //
    arb_clear(lambda);
    arb_clear(logQ);

    compute_zeta_sum(compute_c);
    return 0;
}

int main(int argc, char* argv[])
{
    FILE *file, *outputFile;
    long numbers[MAX_NUMBERS];
    int count = 0;

    file = fopen("input/input.txt", "r");
    if (file == NULL) {
        perror("Error opening file");
        return 1;
    }

    while (fscanf(file, "%ld", &numbers[count]) == 1) {
        count++;
        if (count >= MAX_NUMBERS) {
            printf("Reached maximum number of longs (%d).\n", MAX_NUMBERS);
            break;
        }
    }

    // Get current date and time (without seconds) for filename
    time_t now = time(NULL);
    struct tm *t = localtime(&now);
    char filename[64];
    strftime(filename, sizeof(filename), "output_%Y%m%d%H%M.txt", t);

    // Open output file with date in name
    outputFile = fopen(filename, "w");
    if (outputFile == NULL) {
        perror("Error opening output file");
        return 1;
    }
    printf("started");
    compute_config config;
    init_variables(&config);

    // Write numbers to output file
    for (int i = 0; i < count; i++) {
        if (i% 100 == 0)
        {
            fflush(outputFile);
            printf("%d\n", i);
        }
        long result = compute(&config, numbers[i]);
        fprintf(outputFile, "%ld,%ld\n", numbers[i],result);
    }
    fclose(outputFile);
}
