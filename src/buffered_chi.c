//
// Created by kris on 8/1/25.
//

#include "buffered_chi.h"
#include <stdio.h>
#include <stdlib.h>
#include <flint/ulong_extras.h>

// reads precomputed values of Kronecker symbol
int chi_init(buffered_chi* chi_t, long rows, long cols, char* filename)
{
    chi_t->chi_table = (int*)malloc(rows * cols * sizeof(int));
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
            chi_t->chi_table[i * cols + j] = 1;
        }
        else if (ch == 'N')
        {
            //N signifies nonresidue, i.e. this value is -1
            chi_t->chi_table[i * cols + j] = -1;
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


int chi_val(buffered_chi* chi_t, const long q, const long prime, const long primeIndex)
{
    // compute Kronecker symbol
    int chi = 0;
    if (primeIndex >= chi_t->rows)
    {
        //if we have not precomputed this Kronecker symbol yet, use FLINT
        chi = n_jacobi(q, prime);
    }
    else if (prime == 2 && q % 2 != 0)
    {
        if (q % 8 == 1 || q % 8 == -7)
        {
            chi = 1;
        }
        else
        {
            chi = -1;
        }
    }
    else
    {
        //otherwise use the chi_values array to compute
        long remainder = q % prime;
        if (remainder < 0)
        {
            remainder += prime;
        }
        if (remainder > 0)
        {
            chi = chi_t->chi_table[(primeIndex - 1) * chi_t->cols + remainder - 1];
        }
    }
    return chi;
}
