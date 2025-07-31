#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <flint/arb.h>
#include <flint/ulong_extras.h>
#include <time.h>
#include <unistd.h>

//total number of primes precomputed
#define lenPrime 500000
//FLINT precision to use
#define prec 50
//max number of primes to add to the sum before truncating
#define primeBd 10000
//maximum modulus used
#define qMax 100000000
//dimensions of the Kronecker symbol array
#define rows 10000
#define cols 104729

//reads a precomputed list of primes
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
    printf("Computing Kronecker symbol...\n");
    while(ch != EOF){      
        if (ch == 'R') { //R signifies a residue, i.e. this value is 1
            array[i*cols+j] = 1;
        } else if (ch == 'N') { //N signifies nonresidue, i.e. this value is -1
            array[i*cols +j] = -1;
        } else if (ch == '\n'){
            i++;
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

// loops through all the moduli assigned to one core
void loopQ(long begin, int step, int rank, clock_t start, arb_t lambda, arb_t phi, arb_t O1, int* chi_values, long* primes, FILE* outfile){
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
    arb_t top;
    arb_t bottom;
    arb_t rhs_term_2;
    arb_t div78; //equal to 7/8
    arb_init(div78);
    arb_set_str(div78, "0.875", prec);

    // sets sigma, r
    arb_init(sigma);
    arb_init(logQ);
    arb_init(temp1);
    arb_init(r);
    arb_log_ui(logQ, 10000000000, prec);
    arb_div(r, lambda, logQ, prec);
    arb_add(sigma, one, r, prec);

    for (long q = begin; q <= qMax; q+=step)
    {
        if ((q-begin)/step %100000==0) //print every time the first core does 100000 steps
        {
            printf("%d: Elapsed time when q=%ld: %.3f seconds\n", rank, q, (double)(clock() - start) / CLOCKS_PER_SEC);
        }

        if (q == 0) continue;

        // Calculating log(q)
        arb_init(logq);
        arb_log_ui(logq, abs(q), prec);

        //calculate the partial sum
        arb_init(sum);

        //calculate rhs
        arb_init(rhs);
        arb_init(temp3);
        arb_init(temp4);
        arb_init(temp5);
        arb_init(top);
        arb_init(bottom);
        arb_init(rhs_term_2);

        arb_div(temp3, c, r, prec);
        arb_mul(temp4, r, logq, prec);
        arb_add(temp4, temp4, c, prec);
        arb_div(rhs, temp3, temp4, prec);

        arb_mul(temp5, phi, logq, prec);
        arb_add(rhs, rhs, temp5, prec);
        arb_add(rhs, rhs, O1, prec);

        arb_div(top, c, logq, prec);
        arb_add(top, r, top, prec);
        arb_add(bottom, r, div78, prec);
        arb_mul(bottom, bottom, bottom, prec);
        arb_div(rhs_term_2, top, bottom, prec);
        arb_add(rhs, rhs, rhs_term_2, prec);

        // loop over primes until we exceed primeBd or the inequality is violated
        long prime = primes[0];
        int primeIndex = 0;
        while (primeIndex < primeBd)
        {
            // compute Kronecker symbol
            int chi = 0;
            if(primeIndex>=rows){ //if we have not precomputed this Kronecker symbol yet, use FLINT
                chi = n_jacobi(q, prime);
            }
            else if(prime==2 && q%2 != 0){
                if(q%8==1 || q%8==-7){
                    chi = 1;
                }
                else{
                    chi = -1;
                }
            }
            else { //otherwise use the chi_values array to compute
                long remainder = q%prime;
                if(remainder<0){
                    remainder += prime;
                }
                if(remainder>0){
                    chi = chi_values[(primeIndex-1)*cols+remainder-1];
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

int main(int argc, char** argv)
{

    arb_t lambda;
    arb_t phi;
    arb_t O1; //the O(1)+O(loglog(q)) term
    arb_init(lambda);
    arb_init(phi);
    arb_init(O1);

    //presets:

    //small X (pi(X)~7000)
    arb_set_str(lambda, "1.45", prec);
    arb_set_str(phi, "0.228774", prec);
    arb_set_str(O1, "1.4894", prec);

    //large X (pi(X)~200000)
    // arb_set_str(lambda, "1.3", prec);
    // arb_set_str(phi, "0.23083", prec);
    // arb_set_str(O1, "1.50458", prec);


    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    printf("rank: %d, size: %d\n", rank, size);

    clock_t start, end;
    start = clock(); // start timer on rank 0 only

    long* primes = (long*)malloc(lenPrime * sizeof(long));
    int* chi_values = (int*)malloc(rows * cols * sizeof(int));
    int result = read_kronecker(chi_values); //chi_values[(i-1)*cols+j-1] = z_kronecker(j, primes[i]) when 0<j<primes[i] and primes[i] is odd
    if (read_primes(primes) == 1 || result!=0)
    {
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }

    if (rank == 0)
    {
        printf("Elapsed time after file reading: %.3f seconds\n", (double)(clock() - start) / CLOCKS_PER_SEC);
    }

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
    long begin = 0; //the value q starts at for this core
    int step = 0; //how much q increases by every iteration
    int numEven = (size+1)/2; //number of even-indexed cores
    int numOdd = size/2; //number of odd-indexed cores
    if(rank%2==0){ // if the rank is even, loop through q that are 1 mod 4
        begin = -qMax+1 + 4*(rank/2);
        step = 4*numEven; //evenly distribute the q's that are 1mod4 among the even-indexed cores
    }
    else{ // if the rank is odd, loop through q that are 0 mod 4
        begin = -qMax + 4*(rank/2);
        step = 4*numOdd; //evenly distribute the q's that are 0mod4 among the odd-indexed cores
    }

    loopQ(begin, step, rank, start, lambda, phi, O1, chi_values, primes, outfile);
    if (size==1) { // if there are not multiple threads then we need to do all the calculations on one thread
        loopQ(-qMax, 4, rank, start, lambda, phi, O1, chi_values, primes, outfile);
    }

    fclose(outfile);

    MPI_Finalize();

    if (rank == 0)
    {
        end = clock();
        printf("Elapsed time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    }

    return 0;
}