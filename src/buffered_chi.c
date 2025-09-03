//
// Created by kris on 8/1/25.
//

#include "buffered_chi.h"
#include <stdio.h>
#include <stdlib.h>
#include <flint/ulong_extras.h>

// reads precomputed values of the Kronecker symbol from a file and stores them in chi_t->chi_table
int chi_init(buffered_chi* chi_t, long rows, long cols, char* filename)
{
    chi_t->chi_table = (bool*)calloc(rows * cols, sizeof(bool));
    chi_t->rows = rows;
    chi_t->cols = cols;

    FILE* file = fopen(filename, "rb");
    if (!file)
    {
        perror("Failed to open file");
        return 1;
    }
    int ch = fgetc(file);
    int i = 0;
    int j = 0;
    printf("Computing Kronecker symbol...\n");
    while (ch != EOF)
    {
        if (ch == 'R')
        {
            //R signifies a residue, i.e. this value is 1
            chi_t->chi_table[i * cols + j] = true;
        }
        else if (ch == 'N')
        {
            //N signifies nonresidue, i.e. this value is -1
            chi_t->chi_table[i * cols + j] = false;
        }
        else if (ch == '\n')
        {
            i++;
            j = -1;
        }
        else
        {
            fprintf(stderr, "Invalid character '%c' at line %d, column %d\n", ch, i, j);
            fclose(file);
            return 3;
        }
        ch = fgetc(file);
        j++;
    }

    fclose(file);
    return 0;
}

// calculates the Kronecker symbol (d/prime)_K
int chi_val(buffered_chi* chi_t, const long d, const long prime, const long primeIndex)
{
    // compute Kronecker symbol
    int chi = 0;
    //if we have not precomputed this Kronecker symbol yet, use the FLINT function n_jacobi to calculate
    if (primeIndex >= chi_t->rows)
    {
        chi = n_jacobi(d, prime);
    }
    // the case when prime=2 is not the Legendre symbol
    else if (prime == 2 && d % 2 != 0)
    {
        if (d % 8 == 1 || d % 8 == -7)
        {
            chi = 1;
        }
        else
        {
            chi = -1;
        }
    }
    //otherwise use the chi_values array to compute
    else
    {
        long remainder = d % prime; // the Legendre symbol is prime-periodic so we can reduce d mod prime
        if (remainder < 0)
        {
            remainder += prime;
        }
        if (remainder > 0)
        {
            if (chi_t->chi_table[(primeIndex - 1) * chi_t->cols + remainder - 1]) {
                chi = 1;
            }
            else {
                chi = -1;
            }
        }
    }
    return chi;
}
