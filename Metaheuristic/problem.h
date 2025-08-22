//
// Created by Joan on 12/09/2024.
//
// this .h creates the structure problem_t if it is not done before
#ifndef PROBLEM_H
#define PROBLEM_H
#include <climits>

using namespace std;

#define LITTLE_P

typedef struct {
    int num_S;
    int num_W;
    int *capacity;
    double *fixedcost;
    int *goods;
    double **supplycost;
    int **inc;
    int num_Inc;
    int **f_best_c;
    double **c_best_f;
    int **c_best_f_param;
    int **c_best;
    int **c_best_sc;

} problem_t;


problem_t *create_problem(int, int);
void free_problem(problem_t *);

#endif /* !PROBLEM_H */