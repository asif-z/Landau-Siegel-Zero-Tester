//
// Created by kris on 8/1/25.
//

#ifndef BUFFERED_CHI_H
#define BUFFERED_CHI_H

#include <flint/flint.h>     // for slong, ulong

typedef struct buffered_chi {
    int* chi_table;
    long rows;
    long cols;
} buffered_chi;

int chi_init(buffered_chi* chi_t, long rows, long cols, char* filename);

int chi_val(buffered_chi* chi_t, const slong q, const long prime, const ulong primeIndex);

#endif //BUFFERED_CHI_H


