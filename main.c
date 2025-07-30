#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <flint/arb.h>
#include <flint/long_extras.h>
#include <time.h>
#include <unistd.h>

//total number of primes precomputed
#define lenPrime 500000
//FLINT precision to use
#define prec 50
//max number of primes to add to the sum before truncating
#define primeBd 20000
//maximum modulus used
#define qMax 100000000
//dimensions of the Kronecker symbol array
#define rows 10000
#define cols 104729

//read a list of primes
int read_primes(long* primes)
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

// reads precomputed values of Kronecker symbol
int read_kronecker(int* array) {
    FILE* file = fopen("chi.txt", "rb");
    if (!file) {
        perror("Failed to open file");
        return 1;
    }

    int ch = fgetc(file);
    int i=0;
    int j=0;
    while(ch != EOF){      
        if (ch == 'R') {
            array[i*cols+j] = 1;
        } else if (ch == 'N') {
            array[i*cols +j] = -1;
        } else if (ch == '\n'){
            i++;
            if(i%5000==0){
                printf("Computed Kronecker symbol for %d primes\n", i);
            }
            j = -1;
        }
        else{
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

int main(int argc, char** argv)
{

    arb_t lambda;
    arb_t phi;
    arb_t O1; //the O(1)+O(loglog(q)) term
    arb_init(lambda);
    arb_init(phi);
    arb_init(O1);

    //presets:

    //small X (X~2000)
    arb_set_str(lambda, "2.14", prec);
    arb_set_str(phi, "0.219697", prec);
    arb_set_str(O1, "1.33613", prec);

    //large X (X~50000)
    // arb_set_str(lambda, "1.6", prec);
    // arb_set_str(phi, "0.22675", prec);
    // arb_set_str(O1, "1.38794", prec);


    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    clock_t start, end;
    if (rank == 0)
    {
        start = clock(); // start timer on rank 0 only
    }

    long* primes = (long*)malloc(lenPrime * sizeof(long));
    int* chi_values = malloc(rows * cols * sizeof(int));
    int result = read_kronecker(chi_values); //chi_values[i*cols+j-1] = z_kronecker(j, primes[i]) when 0<j<primes[i]

    if (read_primes(primes) == 1 || result!=0)
    {
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }

    if (rank == 0)
    {
        printf("Elapsed time after file reading: %.3f seconds\n", (double)(clock() - start) / CLOCKS_PER_SEC);
    }

    arb_t c; // zero-free region constant
    arb_init(c);
    arb_set_str(c, "0.1", prec);

    // setting up variables
    arb_t logq;
    arb_t sigma;
    arb_t sum;
    arb_t logp;
    arb_t logQ;
    arb_t p;
    arb_t psigma;
    arb_t zeta_term;
    arb_t one; // equal to 1
    arb_init(one);
    arb_set_ui(one, 1);
    arb_t temp1, temp2, temp3, temp4, temp5; //temp variables for calculations
    arb_t l_term;
    arb_t term;
    arb_t rhs;
    arb_t r;

    // sets sigma, r
    arb_init(sigma);
    arb_init(logQ);
    arb_init(temp1);
    arb_init(r);
    arb_log_ui(logQ, 10000000000, prec);
    arb_div(r, lambda, logQ, prec);
    arb_add(sigma, one, r, prec);

    // Open separate output file per rank to avoid clashes
    char filename[100];
    snprintf(filename, sizeof(filename), "output_rank_%d.csv", rank);
    FILE* outfile = fopen(filename, "w");
    if (outfile == NULL)
    {
        perror("Error opening file");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    //loop through all q
    for (long q = -qMax; q <= qMax; q++)
    {
        if (q % 100000 == 0 & rank == 0)
        {
            printf("Elapsed time when q=%ld: %.3f seconds\n", q, (double)(clock() - start) / CLOCKS_PER_SEC);
        }

        if ((q + qMax) % size != rank) continue; // skip q not assigned to this rank

        if (q == 0) continue;

        // Calculating log(q)
        arb_init(logq);
        arb_log_ui(logq, abs(q), prec);

        if (q % 4 == 0 || q % 4 == 1 || q % 4 == -3)
        {
            //calculate the partial sum
            arb_init(sum);

            //calculate rhs
            arb_init(rhs);
            arb_init(temp3);
            arb_init(temp4);
            arb_init(temp5);

            arb_div(temp3, c, r, prec);
            arb_mul(temp4, r, logq, prec);
            arb_add(temp4, temp4, c, prec);
            arb_div(rhs, temp3, temp4, prec);
            arb_mul(temp5, phi, logq, prec);
            arb_add(rhs, rhs, temp5, prec);
            arb_add(rhs, rhs, O1, prec);

            // loop over primes until we exceed primeBd or the inequality is violated
            long prime = primes[0];
            int primeIndex = 0;
            while (primeIndex < primeBd)
            {
                // compute Kronecker symbol
                int chi = 0;
                if(primeIndex>=rows){ //if we have not precomputed this Kronecker symbol yet, use FLINT
                    chi = z_kronecker(q, prime);
                }
                else { //otherwise use the chi_values array to compute
                    long remainder = q%prime;
                    if(remainder!=0){
                        chi = chi_values[primeIndex*cols+q-1];
                    }
                }
                
                // Calculating log(p)
                arb_init(logp);
                arb_log_ui(logp, prime, prec);

                // Calculating p^sigma
                arb_init(p);
                arb_set_ui(p, prime);
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

                prime = primes[primeIndex++];

                if (primeIndex % 50 == 0 && arb_gt(sum, rhs) == 1)
                {
                    break;
                }
            }

            if (arb_gt(sum, rhs) == 1)
            {
                fprintf(outfile, "%ld,P\n", q);
            }
            else
            {
                fprintf(outfile, "%ld,F\n", q);
            }
        }
    }

    fclose(outfile);

    if (rank == 0)
    {
        end = clock();
        printf("Elapsed time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    }

    MPI_Finalize();

    return 0;
}