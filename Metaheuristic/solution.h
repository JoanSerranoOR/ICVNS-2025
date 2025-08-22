//
// Created by Joan on 17/09/2024.
//

#ifndef SOLUTION_H
#define SOLUTION_H
#include <climits>
#include <set>
#include <stdio.h>
#include <setjmp.h>
#include <sstream>
#include <iterator>

using namespace std;

#define LITTLE_P

typedef struct {
    int num_S;
    int num_W;
    double cost;
    int **assignment_f;
    int **assignment_c;
    int *fopen;
    int *f_assigned;
    int *g_assigned;
    int *c_done;
    int n_open;
    int c_incompleted;
    int **facility_client;
    int **discards;
} solution_t;



solution_t *create_solution(int, int);
void free_solution(solution_t *);

#endif //SOLUTION_H
