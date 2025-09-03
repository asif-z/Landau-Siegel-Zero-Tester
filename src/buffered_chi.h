//
// Created by kris on 8/1/25.
//

#include <stdbool.h>

#ifndef BUFFERED_CHI_H
#define BUFFERED_CHI_H

typedef struct buffered_chi {
    bool* chi_table; //array of precomputed values of the Kronecker symbol
    // the dimensions of the array
    long rows;
    long cols;
} buffered_chi;

int chi_init(buffered_chi* chi_t, long rows, long cols, char* filename);

int chi_val(buffered_chi* chi_t, const long d, const long prime, const long primeIndex);

#endif //BUFFERED_CHI_H


