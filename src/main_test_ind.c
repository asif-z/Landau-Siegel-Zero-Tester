#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "compute.h"

#define MAX_NUMBERS 100000  // Change this as needed

int main(int argc, char* argv[])
{
    FILE *file, *outputFile;
    long numbers[MAX_NUMBERS];
    int count = 0;

    file = fopen("input/output.txt", "r");
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
