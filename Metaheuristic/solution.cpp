#include "solution.h"
#include <stdio.h>
#include <string.h>

solution_t *create_solution(int num_S, int num_W)
{

    auto *solution = new solution_t;

    solution->num_S = num_S;
    solution->num_W = num_W;
    solution->assignment_f = new int*[num_S];
    solution->assignment_c = new int*[num_S];
    solution->facility_client = new int*[num_W];
    solution->fopen = new int[num_W];
    solution->f_assigned = new int[num_W];
    solution->g_assigned = new int[num_S];
    solution->c_done = new int[num_S];
    solution->discards = new int*[num_S];

    for (int i = 0; i < num_S; i++)
    {

        solution->assignment_f[i] = new int[1];
        solution->assignment_c[i] = new int[1];
        solution->discards[i] = new int[2];
    }

   for (int i = 0; i < num_W; i++) {
        solution->facility_client[i] = new int[1];
    }

    return solution;
}


void free_solution(solution_t *solution)
{
    if (solution != NULL) {
        delete []solution;
    }
}

