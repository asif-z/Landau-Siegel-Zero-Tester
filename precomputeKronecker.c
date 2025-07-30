#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <flint/arb.h>
#include <flint/long_extras.h>
#include <time.h>
#include<unistd.h>
#include <stdbool.h>

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
    
    //len of prime
    const long lenPrime = 10000;

    FILE* outfile;
    int count = 0;
    outfile = fopen("chi.txt", "wb");
    if (outfile == NULL)
    {
        perror("Error opening file");
        return 1;
    }

    long* primes = (long*)malloc(lenPrime * sizeof(long));
    read_primes(lenPrime, primes);

    printf("%d\n", primes[lenPrime-1]);
    for (int i=0; i<lenPrime; i++)
    {
        if(i%1000==999){
            printf("done %d\n", i);
        }
        long p = primes[i];
        for(int q=1; q<p; q++){
            int chi = z_kronecker(q, p);
            if(chi==1){
                fprintf(outfile, "R");
            }
            else{
                fprintf(outfile, "N");
            }
        }
        fprintf(outfile, "\n");
    }
    fclose(outfile);

    return 0;
}
