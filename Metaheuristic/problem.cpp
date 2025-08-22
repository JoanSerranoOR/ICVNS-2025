#include <stdio.h>
#include "problem.h"

problem_t *create_problem(int num_S, int num_W)
{

    auto *problem = new problem_t;

    problem->num_S = num_S;
    problem->num_W = num_W;
    problem->capacity = new int[num_W];
    problem->fixedcost = new double[num_W];
    problem->goods = new int[num_S];
    problem->supplycost = new double*[num_S];
    problem->inc = new int*[num_S];
    problem->f_best_c = new int*[num_W];
    problem->c_best_f = new double*[num_S];
    problem->c_best_f_param = new int*[num_S];
    problem->c_best = new int*[num_S];
    problem->c_best_sc = new int*[num_S];


    for (int i = 0; i < num_S; i++)
    {
        problem->inc[i] = new int[num_S];
        problem->supplycost[i] = new double[num_W];
        problem->c_best_f[i] = new double[3];
        problem->c_best_f_param[i] = new int[2];
        problem->c_best[i] = new int[num_W];
        problem->c_best_sc[i] = new int[num_W];
    }
    for (int i = 0; i < num_W; i++)
    {
        problem->f_best_c[i] = new int[num_S];
    }


    return problem;
}


void free_problem(problem_t *problem)
{
    if (problem != NULL) {
        delete[] problem;
    }
}

