# Parallel Processing Assignment 2 - EΠΛ325

## Overview

This repository contains the solutions for the second assignment of the Parallel Processing course (EΠΛ325) during Spring 2023. The assignment focuses on the implementation and analysis of parallel algorithms using OpenMP for an N-body simulation.

## Assignment Description

The assignment consists of implementing three different parallelization methods for an N-body simulation, using OpenMP:

1. **n_body_omp_static**: Implements a static workload distribution among threads.
2. **n_body_omp_dynamic**: Implements a dynamic workload distribution among threads.
3. **n_body_omp_guided**: Implements a guided workload distribution among threads.

The assignment also requires determining the optimal number of threads for each method, analyzing the expected and actual speedup, and calculating the efficiency of each method on the UCY HPC machines.

## Files Included

- `n_body_omp_static.c`: Implementation of the static parallelization method.
- `n_body_omp_dynamic.c`: Implementation of the dynamic parallelization method.
- `n_body_omp_guided.c`: Implementation of the guided parallelization method.
- `123456_HW2.pdf`: Detailed report including analysis, results, and observations.
- `123456_HW2.xlsx`: Spreadsheet with the data and graphs generated for the report.
- `README.md`: This file, providing an overview of the assignment and repository.
- Other supporting files used in the implementation and testing of the solutions.

## Compilation Instructions

To compile the code with OpenMP support, use the following command:

```bash
gcc -Wall -Werror -fopenmp -O3 -o n_body_static n_body_omp_static.c
gcc -Wall -Werror -fopenmp -O3 -o n_body_dynamic n_body_omp_dynamic.c
gcc -Wall -Werror -fopenmp -O3 -o n_body_guided n_body_omp_guided.c
