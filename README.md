# Landau-Siegel Zero Tester

This is part of an undergraduate research project at University of Toronto, Math Department that 
aims to verify the non-existence of Landau-Siegel zero in the region $\beta_1 > 1-\frac{c}{\log(q)}$
for primitive quadratic characters mod $q$, with $0<q<q_{max}$. This project was funded by NSERC and 
UTEA Undergraduate Student Research Awards. This research is also supported by Compute Canada 
([alliancecan.ca](https://alliancecan.ca/en)), where most of our computation was performed.

## Table of Contents
1. [background](#background)
2. [Requirements](#requirements)
3. [Usage](#usage)
4. [Directory Structure](#directory-structure)

## Background

Let $\chi$ be a Dirichlet character mod $q$,
$$L(s,\chi) := \sum_{n=1}^\infty \frac{\chi(n)}{n^s}, \quad \Re(s) > 0,$$
be the associated Dirichlet $L$-function, and $s = \sigma + it$ be a complex number.

Landau proved that there is a constant $c > 0$ such that for all $q$,
$\prod_\chi L(s,\chi)$ has at most one zero in the region

$$
\sigma > 1 - \frac{c}{\log(q(|t|+4))}
$$

where $\chi$ ranges over all Dirichlet characters of modulus $q$.  
If such a zero exists, then it is necessarily real and the associated character $\chi$ is quadratic.

Such a zero, if it exists, is called a **Landau–Siegel zero** or an **exceptional zero**.  

Running this program on Compute Canada's Nibi clusters, we computationally verified that Landau–Siegel zeros do not 
exist for any primitive quadratic character of modulus $q \le 10^{9}$. In particular, we proved the 
following theorem.

**Theorem 1.**  
Let $\chi$ be a primitive quadratic character of modulus $q \leq 10^{9}$. Then $L(s,\chi)$ has no real zeros $\beta_1$ in the region

$$
\beta_1 > 1 - \frac{c}{\log q} 
$$

with $c = \tfrac{1}{10}$.

## Requirements

To build and run the computational verification code, you will need the following system dependencies:

- [**FLINT**](http://flintlib.org/) (Fast Library for Number Theory) ≥ 2.9
- **CMake** ≥ 3.20 (for building the project)
- **MPI** (Message Passing Interface, e.g. OpenMPI or MPICH) for parallel execution

Make sure these are installed and available in your system path.

## Usage

After installing the requirements and building the project, you can run the computations as follows.

### 1. Precompute chi

Run 'misc/precompute_kronecker.c' to generate `chi'txt`. After the file is generated, move it to `/input` dirctory.

### 2. Choose the version

Edit `src/CMakeLists.txt` to choose the version of main program to run.

### 3. Build the project
From the repository root:  
```bash
mkdir build
cd build
cmake ..
cmake --build .
```

### 4. Run with MPI
From the repository root use the following mpi command to run, change 10 to desired number of MPI processes.
```bash
mpirun -np 10 ./build/src/mpiTest
```

### 5. Using `main_by_file.c`
If one wants to run the verification on a list of $q$ in a file, name the file to `input.txt` and store it in `/input` 
dirctory.

## Directory Structure

All source code is stored in the directory `/src`. All input files are stored in `/input` dirctory. 
`/test` dirctory contains all unit tests for this project and `/misc` conatins all other helping code and scripts.
