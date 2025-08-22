#include "problem.h"
#include "solution.h"
#include "heuristic.h"
#include "globals.h"
//#include "plotwidget.h"
//#include "plotwidget_detail.h"
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include <chrono>
#include <cstring>
#include <fstream>
#include <io.h>
#include <locale>
#include <cstdlib>
/*#include <qmath.h>
#include <QGridLayout>
#include <QDebug>
#include <QMainWindow>
#include <QtCharts/QChartView>
#include <QApplication>
#include <QWidget>
#include <QVBoxLayout>
#include <utility>*/
using namespace std;

void guardarSolucion3b(solution_t *solucion, const string &file_name, int S, int W);
Heuristic::Heuristic(problem_t *problem)
{
    this->problem = problem;
}

Heuristic::~Heuristic()
{
}

void Heuristic::initialize() {
}

void Heuristic::init_solution(const problem_t *problem, solution_t *ini_sol) {
    int S = problem->num_S;
    int W = problem->num_W;
    ini_sol->num_S = S;
    ini_sol->num_W = W;

    for (int j = 0; j < W; j++)
    {
        ini_sol->fopen[j] = 0;
        ini_sol->f_assigned[j] = 0;
        ini_sol->facility_client[j][0] = 0;
    }

    for (int j = 0; j < S; j++)
    {
        ini_sol->c_done[j] = 0;
        ini_sol->g_assigned[j] = 0;

    }

    for (int i = 0; i < S; i++) {
        ini_sol->assignment_f[i][0] = 0;
        ini_sol->assignment_c[i][0] = 0;
    }


    ini_sol->cost =0;
    ini_sol->n_open = 0;
    ini_sol->c_incompleted = 0;

    for (int i = 0; i < S; i++) {
        ini_sol->discards[i][0] = 0;
        ini_sol->discards[i][1] = 0;
    }

}

int minimo(int min, int a) {
    if (min > a) {
        min = a;
    }
    return min;
}

int maximo(int max, int a) {
    if (max < a) {
        max = a;
    }
    return max;
}

double minimo_double(double min, double a) {
    if (min > a) {
        min = a;
    }
    return min;
}

double maximo_double(double max, double a) {
    if (max < a) {
        max = a;
    }
    return max;
}

vector<vector<int>> Heuristic::buildRCL1(const problem_t *problem, solution_t *sol) {
    std::vector<std::vector<double>> RCL1;
    std::vector<std::vector<double>> RCL;
    std::vector<std::vector<int>> RCL_return;
    double a;
    double min = 1000000000.0;
    double max = 0.0;

    for (int i = 0; i < problem->num_S; i++) {
        int cliente = i;


        if ((sol->c_done[cliente] < 2)) {


            if (sol->discards[cliente][0] == 0) {

                int iter = 0;
                while (iter < problem->num_W) {
                    int facility = problem->c_best[cliente][iter];


                    if (sol->f_assigned[facility] == problem->capacity[facility]) {
                        iter++;
                    }else {

                        if (c_best_mode == 1) {
                            a = (problem->goods[cliente]-sol->g_assigned[cliente])*problem->supplycost[cliente][facility];
                        }else if (c_best_mode == 2) {
                            a = (problem->goods[cliente]-sol->g_assigned[cliente])*problem->supplycost[cliente][facility] + problem->fixedcost[facility];
                        }else {
                            a = (problem->goods[cliente]-sol->g_assigned[cliente])*problem->supplycost[cliente][facility] + problem->fixedcost[facility]/problem->capacity[facility];
                        }

                        RCL1.push_back({(double)cliente,(double)facility,a});
                        break;
                    }
                }

            }else {


                int l = sol->discards[cliente][1];

                for (int j = 0; j < problem->num_W; j++) {
                    int suma = 0;
                    for (int k = 0; k < l; k++) {
                        if (problem->c_best[cliente][j] == sol->discards[cliente][k+2]) {

                            suma++;
                            break;
                        }
                    }
                    if (suma == 0) {
                        int facility = problem->c_best[cliente][j];

                        if (sol->f_assigned[facility] < problem->capacity[facility]) {

                            if (c_best_mode == 1) {
                                a = (problem->goods[cliente]-sol->g_assigned[cliente])*problem->supplycost[cliente][facility];
                            }else if (c_best_mode == 2) {
                                a = (problem->goods[cliente]-sol->g_assigned[cliente])*problem->supplycost[cliente][facility] + problem->fixedcost[facility];
                            }else {
                                a = (problem->goods[cliente]-sol->g_assigned[cliente])*problem->supplycost[cliente][facility] + problem->fixedcost[facility]/problem->capacity[facility];
                            }
                            RCL1.push_back({(double)cliente,(double)facility,a});
                            break;
                        }

                    }
                }
            }


            min = minimo_double(min, a);
            max = maximo_double(max, a);
        }
    }

    for (int j = 0; j < RCL1.size(); j++) {

        double cliente = RCL1[j][0];
        double facility = RCL1[j][1];
        double coste = RCL1[j][2];
        if ((coste >=  max - alpha*(max-min))) {
            RCL.push_back({cliente,facility,coste});
        }
    }
    std::sort(RCL.begin(), RCL.end(), [](const std::vector<double>& a, const std::vector<double>& b) {
        return a[2] > b[2];
        });

    int msize = RCL.size();
    for (int i = 0; i < msize; i++) {
        RCL_return.push_back({(int)RCL[i][0],(int)RCL[i][1]});
    }

    return RCL_return;
}

vector<vector<int>> Heuristic::CostPerCustomer(const problem_t *problem, const solution_t *sol) {
    std::vector<vector<int>> RCL;
    int facility;
    int coste;
    for (int i = 0; i < problem->num_S; i++) {
        coste = 0;
        int cliente = i;

        for (int j = 1; j < sol->assignment_f[cliente][0]+1; j++) {
            facility = sol->assignment_f[cliente][j];
            coste = coste + sol->assignment_c[cliente][j]*problem->supplycost[cliente][facility];
        }

        RCL.push_back({cliente,coste});
    }

    std::sort(RCL.begin(), RCL.end(), [](const std::vector<int>& a, const std::vector<int>& b) {
        return a[1] > b[1];
        });

    return RCL;
}

vector<vector<int>> Heuristic::CostPerCustomerInc(const problem_t *problem, const solution_t *sol) {
    std::vector<vector<int>> RCL;
    int facility;
    int coste;
    for (int i = 0; i < problem->num_S; i++) {
        coste = 0;
        int cliente = i;

        int cond = 0;
        for (int r = 0; r < problem->num_S; r++) {
            if (cliente != r) {
                if (problem->inc[cliente][r] == 1) {
                    cond++;
                    break;
                }
            }
        }

        if (cond > 0) {
            for (int j = 1; j < sol->assignment_f[cliente][0]+1; j++) {
                facility = sol->assignment_f[cliente][j];
                coste = coste + sol->assignment_c[cliente][j]*problem->supplycost[cliente][facility];
            }

            RCL.push_back({cliente,coste});
        }

    }

    std::sort(RCL.begin(), RCL.end(), [](const std::vector<int>& a, const std::vector<int>& b) {
        return a[1] > b[1];
        });

    return RCL;
}

void Heuristic::DropCustomer(const problem_t *problem, solution_t *solution, int cliente, int facility_old) {
    int pos_facility;
    for (int i = 1; i < solution->assignment_f[cliente][0]+1; i++) {
        if (solution->assignment_f[cliente][i] == facility_old) {
            pos_facility = i;
            break;
        }
    }

    int q = solution->assignment_c[cliente][pos_facility];

    solution->g_assigned[cliente] = solution->g_assigned[cliente] - q;
    solution->f_assigned[facility_old] = solution->f_assigned[facility_old] - q;

    if (solution->c_done[cliente] == 2 && q > 0) {
        solution->c_incompleted--;
    }

    if (solution->g_assigned[cliente] == 0) {
        solution->c_done[cliente] = 0;
    }else {
        solution->c_done[cliente] = 1;
    }


    solution->cost = solution->cost - q*problem->supplycost[cliente][facility_old];
    int open = solution->fopen[facility_old];
    if (open == 1 && solution->f_assigned[facility_old] == 0) {
        solution->fopen[facility_old] = 0;
        solution->n_open--;
        solution->cost = solution->cost - problem->fixedcost[facility_old];
    }

    int l = solution->assignment_f[cliente][0];
    vector<int> assignment_c;
    vector<int> assignment_f;

    assignment_c.push_back(l-1);
    assignment_f.push_back(l-1);

    for (int i = 1; i < l+1; i++) {
        if (solution->assignment_f[cliente][i] != facility_old) {
            assignment_c.push_back(solution->assignment_c[cliente][i]);
            assignment_f.push_back(solution->assignment_f[cliente][i]);
        }
    }
    delete [] solution->assignment_c[cliente];
    solution->assignment_c[cliente] = new int [l];
    delete [] solution->assignment_f[cliente];
    solution->assignment_f[cliente] = new int [l];

    for (int i = 0; i < l; i++) {
        solution->assignment_c[cliente][i] = assignment_c[i];
        solution->assignment_f[cliente][i] = assignment_f[i];
    }

    int l2 = solution->facility_client[facility_old][0];
    vector <int> v_facilities;

    v_facilities.push_back(l2-1);
    for (int i = 1; i < l2+1; i++) {
        if (solution->facility_client[facility_old][i] != cliente) {
            int r = solution->facility_client[facility_old][i];
            v_facilities.push_back(r);
        }
    }
    delete [] solution->facility_client[facility_old];
    solution->facility_client[facility_old] = new int [l2];

    for (int j = 0; j < l2; j++) {
        solution->facility_client[facility_old][j] = v_facilities[j];
    }
}

void Heuristic::DropCustomer_Split(const problem_t *problem, solution_t *solution, int cliente, int facility_old, int q, int splits) {
    if (q > 0) {
        int pos_facility;

    for (int i = 1; i < solution->assignment_f[cliente][0]+1; i++) {
        if (solution->assignment_f[cliente][i] == facility_old) {
            pos_facility = i;
            break;
        }
    }


    solution->g_assigned[cliente] = solution->g_assigned[cliente] - q;
    solution->f_assigned[facility_old] = solution->f_assigned[facility_old] - q;


    if (solution->c_done[cliente] == 2 && q > 0) {
        solution->c_incompleted--;
    }

    if (solution->g_assigned[cliente] == 0) {
        solution->c_done[cliente] = 0;
    }else {
        solution->c_done[cliente] = 1;
    }


    solution->cost = solution->cost - q*problem->supplycost[cliente][facility_old];
    int open = solution->fopen[facility_old];
    if (open == 1 && solution->f_assigned[facility_old] == 0) {
        solution->fopen[facility_old] = 0;
        solution->n_open--;
        solution->cost = solution->cost - problem->fixedcost[facility_old];
    }

    if (splits == 1) {

        int l = solution->assignment_f[cliente][0];

        vector<int> assignment_c;
        vector<int> assignment_f;

        assignment_c.push_back(l-1);
        assignment_f.push_back(l-1);

        for (int i = 1; i < l+1; i++) {
            if (solution->assignment_f[cliente][i] != facility_old) {
                assignment_c.push_back(solution->assignment_c[cliente][i]);
                assignment_f.push_back(solution->assignment_f[cliente][i]);
            }
        }

        delete [] solution->assignment_c[cliente];
        solution->assignment_c[cliente] = new int [l];
        delete [] solution->assignment_f[cliente];
        solution->assignment_f[cliente] = new int [l];

        for (int i = 0; i < l; i++) {
            solution->assignment_c[cliente][i] = assignment_c[i];
            solution->assignment_f[cliente][i] = assignment_f[i];
        }

        int l2 = solution->facility_client[facility_old][0];
        vector <int> v_facilities;

        v_facilities.push_back(l2-1);
        for (int i = 1; i < l2+1; i++) {
            if (solution->facility_client[facility_old][i] != cliente) {
                int r = solution->facility_client[facility_old][i];
                v_facilities.push_back(r);
            }
        }
        delete [] solution->facility_client[facility_old];
        solution->facility_client[facility_old] = new int [l2];

        for (int j = 0; j < l2; j++) {
            solution->facility_client[facility_old][j] = v_facilities[j];
        }
    }else {
        solution->assignment_c[cliente][pos_facility] = solution->assignment_c[cliente][pos_facility] - q;
    }
    }

}

void Heuristic::InsertCustomer_2(const problem_t *problem, solution_t *sol, int cliente, int facility, int done) {
    int l, pos2;
    int f_open = sol->fopen[facility];

    int suma_valor = minimo(problem->goods[cliente]-sol->g_assigned[cliente],problem->capacity[facility]-sol->f_assigned[facility]);

    if ((sol->f_assigned[facility] + suma_valor <= problem->capacity[facility]) && (done < 2) && (suma_valor != 0)) {
    sol->cost = sol->cost + suma_valor*problem->supplycost[cliente][facility];
    sol->fopen[facility] = 1;
    sol->f_assigned[facility] = sol->f_assigned[facility] + suma_valor;
    sol->g_assigned[cliente] = sol->g_assigned[cliente] + suma_valor;
    if (sol->g_assigned[cliente] == problem->goods[cliente]) {
        sol->c_done[cliente] = 2;
        if (done < 2 && sol->c_done[cliente] == 2) {
            sol->c_incompleted++;
        }
    }else if (sol->g_assigned[cliente] > 0) {
        sol->c_done[cliente] = 1;
    }

    if (f_open == 0 && sol->fopen[facility] == 1) {
        sol->cost = sol->cost + problem->fixedcost[facility];
        sol->n_open++;
    }

    l = sol->assignment_f[cliente][0];
    std::vector<int> matrix_f;
    std::vector<int> matrix_c;

    for (int i = 0; i < l+1; i++) {
        matrix_f.push_back(sol->assignment_f[cliente][i]);
        matrix_c.push_back(sol->assignment_c[cliente][i]);
    }
    auto start = matrix_f.begin()+1;
    auto end = start + l;

    auto pos = std::find(start, end, facility);
    if (pos != end) {
        pos2 = distance(start, pos);
        matrix_c[pos2+1] = matrix_c[pos2+1] + suma_valor;

        delete [] sol->assignment_c[cliente];
        sol->assignment_c[cliente] = new int [l+1];

        for (int i = 0; i < l+1; i++) {
            sol->assignment_c[cliente][i] = matrix_c[i];
        }
    }else {
        matrix_f.push_back(facility);
        matrix_c.push_back(suma_valor);

        delete [] sol->assignment_f[cliente];
        delete [] sol->assignment_c[cliente];
        sol->assignment_f[cliente] = new int [l+2];
        sol->assignment_c[cliente] = new int [l+2];

        sol->assignment_f[cliente][0] = matrix_f[0]+1;
        sol->assignment_c[cliente][0] = matrix_c[0]+1;
        for (int i = 1; i < l+2; i++) {
            sol->assignment_f[cliente][i] = matrix_f[i];
            sol->assignment_c[cliente][i] = matrix_c[i];
        }
    }
    int longitud = sol->facility_client[facility][0];
    std::vector<int> v_temporal(longitud+2);

    for (int i = 0; i < longitud+1; i++) {
            v_temporal[i] = sol->facility_client[facility][i];
    }
    v_temporal[longitud+1] = -1;


    if (longitud != 0) {
        auto start = v_temporal.begin()+1;
        auto end = v_temporal.end();
        pos = find(start, end, cliente);
        if (pos == end) {
            v_temporal[0] = v_temporal[0]+1;

            v_temporal[longitud+1] = cliente;

            delete [] sol->facility_client[facility];

            sol->facility_client[facility] = new int [longitud+2];

            for (int i = 0; i < longitud+2; i++) {
                sol->facility_client[facility][i] = v_temporal[i];
            }


        }
    }else {
        v_temporal[0] = v_temporal[0]+1;

        v_temporal[longitud+1] = cliente;

        delete [] sol->facility_client[facility];

        sol->facility_client[facility] = new int [longitud+2];

        for (int i = 0; i < longitud+2; i++) {
            sol->facility_client[facility][i] = v_temporal[i];
        }
    }

    }else if (problem->capacity[facility]==sol->f_assigned[facility]){
        vector<int> v_temp;

        int l = sol->discards[cliente][1];
        for (int i = 0; i < l+2; i++) {
            v_temp.push_back(sol->discards[cliente][i]);
        }
        auto start = v_temp.begin()+1;
        auto end = v_temp.end();

        auto pos = find(start, end, facility);

        if (l == 0 | pos == end) {
            v_temp[1]++;
            v_temp.push_back(facility);
            delete [] sol->discards[cliente];
            sol->discards[cliente] = new int[l+3];

            for (int i = 0; i < l+3; i++) {
                sol->discards[cliente][i] = v_temp[i];
            }
        }

    }
}

void Heuristic::InsertCustomer_2q(const problem_t *problem, solution_t *sol, int cliente, int facility, int suma_valor ,int done) {
    int l, pos2;
    int f_open = sol->fopen[facility];

if ((sol->f_assigned[facility] + suma_valor <= problem->capacity[facility]) && (done < 2) && (suma_valor != 0)) {
    sol->cost = sol->cost + suma_valor*problem->supplycost[cliente][facility];
    sol->fopen[facility] = 1;
    sol->f_assigned[facility] = sol->f_assigned[facility] + suma_valor;
    sol->g_assigned[cliente] = sol->g_assigned[cliente] + suma_valor;

    if (sol->g_assigned[cliente] == problem->goods[cliente]) {
        sol->c_done[cliente] = 2;
        if (done < 2 && sol->c_done[cliente] == 2) {
            sol->c_incompleted++;
        }
    }else if (sol->g_assigned[cliente] > 0) {
        sol->c_done[cliente] = 1;
    }

    if (f_open == 0 && sol->fopen[facility] == 1) {
        sol->cost = sol->cost + problem->fixedcost[facility];
        sol->n_open++;
    }

    l = sol->assignment_f[cliente][0];
    std::vector<int> matrix_f;
    std::vector<int> matrix_c;

    for (int i = 0; i < l+1; i++) {
        matrix_f.push_back(sol->assignment_f[cliente][i]);
        matrix_c.push_back(sol->assignment_c[cliente][i]);
    }
    auto start = matrix_f.begin()+1;
    auto end = start + l;

    auto pos = std::find(start, end, facility);
    if (pos != end) {
        pos2 = distance(start, pos);
        matrix_c[pos2+1] = matrix_c[pos2+1] + suma_valor;

        delete [] sol->assignment_c[cliente];
        sol->assignment_c[cliente] = new int [l+1];

        for (int i = 0; i < l+1; i++) {
            sol->assignment_c[cliente][i] = matrix_c[i];
        }
    }else {
        matrix_f.push_back(facility);
        matrix_c.push_back(suma_valor);

        delete [] sol->assignment_f[cliente];
        delete [] sol->assignment_c[cliente];
        sol->assignment_f[cliente] = new int [l+2];
        sol->assignment_c[cliente] = new int [l+2];

        sol->assignment_f[cliente][0] = matrix_f[0]+1;
        sol->assignment_c[cliente][0] = matrix_c[0]+1;
        for (int i = 1; i < l+2; i++) {
            sol->assignment_f[cliente][i] = matrix_f[i];
            sol->assignment_c[cliente][i] = matrix_c[i];
        }
    }

    int longitud = sol->facility_client[facility][0];
    std::vector<int> v_temporal(longitud+2);

    for (int i = 0; i < longitud+1; i++) {
            v_temporal[i] = sol->facility_client[facility][i];
    }
    v_temporal[longitud+1] = -1;


    if (longitud != 0) {
        auto start = v_temporal.begin()+1;
        auto end = v_temporal.end();
        pos = find(start, end, cliente);
        if (pos == end) {
            v_temporal[0] = v_temporal[0]+1;

            v_temporal[longitud+1] = cliente;

            delete [] sol->facility_client[facility];

            sol->facility_client[facility] = new int [longitud+2];

            for (int i = 0; i < longitud+2; i++) {
                sol->facility_client[facility][i] = v_temporal[i];
            }


        }
    }else {

        v_temporal[0] = v_temporal[0]+1;

        v_temporal[longitud+1] = cliente;

        delete [] sol->facility_client[facility];

        sol->facility_client[facility] = new int [longitud+2];

        for (int i = 0; i < longitud+2; i++) {
            sol->facility_client[facility][i] = v_temporal[i];
        }
    }

    }else if (problem->capacity[facility]==sol->f_assigned[facility]){

        vector<int> v_temp;

        int l = sol->discards[cliente][1];
        for (int i = 0; i < l+2; i++) {
            v_temp.push_back(sol->discards[cliente][i]);
        }
        auto start = v_temp.begin()+1;
        auto end = v_temp.end();

        auto pos = find(start, end, facility);

        if (l == 0 | pos == end) {
            v_temp[1]++;
            v_temp.push_back(facility);
            delete [] sol->discards[cliente];
            sol->discards[cliente] = new int[l+3];

            for (int i = 0; i < l+3; i++) {
                sol->discards[cliente][i] = v_temp[i];
            }
        }

    }

}

void Heuristic::FillFacility_2b(const problem_t *problem, solution_t *sol, int facility) {
    int cliente, done;
    int iter = 0;

    while ((sol->f_assigned[facility] < problem->capacity[facility]) && (iter < problem->num_S)) {
        cliente = problem->f_best_c[facility][iter];

        done = sol->c_done[cliente];

        if (problem->supplycost[cliente][facility]<= Dmax) {
            if (done < 2) {
                int suma_inc = 0;
                int n_max;

                for (int i = 0; i < problem->num_S; i++) {
                    n_max = sol->assignment_f[i][0];
                    for (int j = 0; j < n_max; j++) {
                        if (sol->assignment_f[i][j+1] == facility && i != cliente && (problem->inc[cliente][i] > 0 || problem->inc[i][cliente] > 0)) {
                            suma_inc++;
                            break;
                        }
                    }
                    if (suma_inc > 0) {
                        break;
                    }
                }

                if (suma_inc == 0) {

                    if (problem->capacity[facility] - sol->f_assigned[facility] >= problem->goods[cliente] - sol->g_assigned[cliente]) {
                        InsertCustomer_2(problem, sol, cliente, facility, done);
                    }
                }
                else {

                    vector<int> v_temp;

                    int l = sol->discards[cliente][1];
                    for (int i = 0; i < l+2; i++) {
                        v_temp.push_back(sol->discards[cliente][i]);
                    }
                    auto start = v_temp.begin()+1;
                    auto end = v_temp.end();

                    auto pos = find(start, end, facility);

                    if (l == 0 | pos == end) {
                        v_temp[1]++;
                        v_temp.push_back(facility);
                        delete [] sol->discards[cliente];
                        sol->discards[cliente] = new int[l+3];

                        for (int i = 0; i < l+3; i++) {
                            sol->discards[cliente][i] = v_temp[i];
                        }
                    }

                }
            }
        }else {
            break;
        }

        iter++;
    }
}

void Heuristic::CompleteSolb_shake(const problem_t *problem, solution_t *complete_sol) {
    int iter = 0;
    std::vector<std::vector<int>> RCL;
    int random_number, cliente, facility, f_open, done;

    for (int i = 0 ; i < problem->num_S; i++) {
        delete [] complete_sol->discards[i];
        complete_sol->discards[i] = new int[2];

        complete_sol->discards[i][0] = 0;
        complete_sol->discards[i][1] = 0;

    }


    while (complete_sol->c_incompleted < problem->num_S){

        RCL = buildRCL1(problem, complete_sol);

        int max_row = RCL.size();
        if (max_row > 0) {

            std::uniform_int_distribution<int> dist(0, max_row-1);

            random_number = dist(rng);

            cliente = RCL[random_number][0];
            facility = RCL[random_number][1];

            if (problem->supplycost[cliente][facility] <= Dmax) {
                f_open = complete_sol->fopen[facility];
                done = complete_sol->c_done[cliente];



                int suma_inc = 0;
                int n_max ;
                for (int i = 0; i < problem->num_S; i++) {
                    n_max = complete_sol->assignment_f[i][0];
                    for (int j = 0; j < n_max; j++) {
                        if (complete_sol->assignment_f[i][j+1] == facility && i != cliente && (problem->inc[cliente][i] > 0 || problem->inc[i][cliente] > 0)) {
                            suma_inc++;
                            break;
                        }
                    }
                    if (suma_inc > 0) {
                        break;
                    }
                }

                if (suma_inc == 0 && (problem->capacity[facility]-complete_sol->f_assigned[facility] >=
                    problem->goods[cliente]-complete_sol->g_assigned[cliente])) {
                    InsertCustomer_2(problem, complete_sol, cliente, facility, done);
                    FillFacility_2b(problem, complete_sol, facility);
                }
                else {
                    complete_sol->discards[cliente][0] = 1;

                    vector<int> v_temp;
                    int l = complete_sol->discards[cliente][1];

                    for (int i = 0; i < l+2; i++) {
                        v_temp.push_back(complete_sol->discards[cliente][i]);
                    }
                    v_temp[1]++;
                    v_temp.push_back(facility);
                    delete [] complete_sol->discards[cliente];
                    complete_sol->discards[cliente] = new int[l+3];

                    for (int i = 0; i < l+3; i++) {
                        complete_sol->discards[cliente][i] = v_temp[i];
                    }
                }

            }else {
                complete_sol->discards[cliente][0] = 1;

                vector<int> v_temp;
                int l = complete_sol->discards[cliente][1];

                for (int i = 0; i < l+2; i++) {
                    v_temp.push_back(complete_sol->discards[cliente][i]);
                }
                v_temp[1]++;
                v_temp.push_back(facility);
                delete [] complete_sol->discards[cliente];
                complete_sol->discards[cliente] = new int[l+3];

                for (int i = 0; i < l+3; i++) {
                    complete_sol->discards[cliente][i] = v_temp[i];
                }
            }
        }else {
            break;
        }
        iter++;
    }

}

void Heuristic::CompleteSolb(const problem_t *problem, solution_t *complete_sol) {
    int iter = 0;
    std::vector<std::vector<int>> RCL;
    int random_number, cliente, facility, f_open, done;

    while (complete_sol->c_incompleted < problem->num_S){
        RCL = buildRCL1(problem, complete_sol);

        int max_row = RCL.size();

        if (max_row > 0) {

            std::uniform_int_distribution<int> dist(0, max_row-1);

            random_number = dist(rng);

            cliente = RCL[random_number][0];
            facility = RCL[random_number][1];

            if (problem->supplycost[cliente][facility]<=Dmax) {
                f_open = complete_sol->fopen[facility];
                done = complete_sol->c_done[cliente];

                int suma_inc = 0;
                int n_max ;
                for (int i = 0; i < problem->num_S; i++) {
                    n_max = complete_sol->assignment_f[i][0];
                    for (int j = 0; j < n_max; j++) {
                        if (complete_sol->assignment_f[i][j+1] == facility && i != cliente && (problem->inc[cliente][i] > 0 || problem->inc[i][cliente] > 0)) {
                            suma_inc++;
                            break;
                        }
                    }
                    if (suma_inc > 0) {
                        break;
                    }
                }

                if (suma_inc == 0 && (problem->capacity[facility]-complete_sol->f_assigned[facility] >=
                    problem->goods[cliente]-complete_sol->g_assigned[cliente])) {
                    InsertCustomer_2(problem, complete_sol, cliente, facility, done);
                    FillFacility_2b(problem, complete_sol, facility);
                }
                else {
                    complete_sol->discards[cliente][0] = 1;

                    vector<int> v_temp;
                    int l = complete_sol->discards[cliente][1];

                    for (int i = 0; i < l+2; i++) {
                        v_temp.push_back(complete_sol->discards[cliente][i]);
                    }
                    v_temp[1]++;
                    v_temp.push_back(facility);
                    delete [] complete_sol->discards[cliente];
                    complete_sol->discards[cliente] = new int[l+3];

                    for (int i = 0; i < l+3; i++) {
                        complete_sol->discards[cliente][i] = v_temp[i];
                    }
                }
            }else {
                complete_sol->discards[cliente][0] = 1;

                vector<int> v_temp;
                int l = complete_sol->discards[cliente][1];

                for (int i = 0; i < l+2; i++) {
                    v_temp.push_back(complete_sol->discards[cliente][i]);
                }
                v_temp[1]++;
                v_temp.push_back(facility);
                delete [] complete_sol->discards[cliente];
                complete_sol->discards[cliente] = new int[l+3];

                for (int i = 0; i < l+3; i++) {
                    complete_sol->discards[cliente][i] = v_temp[i];
                }
            }


        }else {
            break;
        }
        iter++;
    }




}

void Heuristic::Constructive_8b(const problem_t *problem, solution_t *ini_sol) {
    init_solution(problem, ini_sol);

    int cont = 0;
    int fac, cli, cliente, facility, done, random_number, fac2, fac3;
    int max_cost;

    vector<int> open_f;
    vector<int> open_f2;
    vector<vector<int> > count_closestfac;

    for (int i = 0; i < problem->num_W; i++) {
        count_closestfac.push_back({i, 0});
    }

    for (int i = 0; i < problem->num_S; i++) {
        fac = problem->c_best[i][0];
        count_closestfac[fac][1]++;
    }
    for (int i = 0; i < problem->num_W; i++) {
        count_closestfac[i][1]++;
    }

    int facopen = 0;

    for (int j = 0; j < problem->num_W; j++) {
        facopen = facopen + ini_sol->fopen[j];
    }

    int fac_to_open = minimo(problem->num_W-facopen, ceil(ini_open* problem->num_W));

    while (cont < fac_to_open) {
        std::vector<int> v_facilities;
        for (const auto &row: count_closestfac) {
            v_facilities.push_back(row[1]);
        }
        std::discrete_distribution<int> distribution(v_facilities.begin(), v_facilities.end());
        int cond = 0;
        while (cond == 0) {
            facility = distribution(rng);

            if (open_f.size() == 0) {
                cond++;

            } else {
                auto it2 = find(open_f.begin(), open_f.end(), facility);
                if ((it2 == open_f.end()) && ini_sol->fopen[facility] == 0) {
                    cond++;
                }
            }
        }


        ini_sol->fopen[facility] = 1;
        ini_sol->n_open++;
        ini_sol->cost = ini_sol->cost + problem->fixedcost[facility];
        open_f.push_back(facility);

        open_f2.push_back(facility);
        cont++;
        count_closestfac[facility][1] = 0;
    }

    cont = 0;
    while ((cont < ini_sol->n_open) && (ini_sol->c_incompleted < problem->num_S)) {
        vector<vector<int>> count_lessDmax;

        for (int i = 0; i < problem->num_S; i++) {
            count_lessDmax.push_back({i, 0});
        }

        for (int i = 0; i < problem->num_S; i++) {
            if (ini_sol->c_done[i] == 2) {
                continue;
            }else {
                int pos = 0;

                while(pos < problem->num_W) {
                    int fac = problem->c_best_sc[i][pos];

                    if (problem->supplycost[i][fac] > Dmax) {
                        break;
                    }
                    pos++;
                }

                count_lessDmax[i][1] = pos+1;


            }
        }

        std::vector<int> v_facilities2;
        for (const auto &row: count_lessDmax) {
            v_facilities2.push_back(row[1]);
        }

        for (int i = 0; i < problem->num_S; i++) {
            if (v_facilities2[i] != 0) {
                v_facilities2[i] = 100.0/v_facilities2[i];
            }

        }

        std::discrete_distribution<int> disc_dist(v_facilities2.begin(), v_facilities2.end());

        int c1 = disc_dist(rng);
        int fac_fill;
        for (int j1 = 0; j1 < problem->num_W; j1++) {

            fac_fill = problem->c_best_sc[c1][j1];

            if ((ini_sol->f_assigned[fac_fill] == 0) && (ini_sol->fopen[fac_fill] == 1)) {
                InsertCustomer_2(problem, ini_sol, c1, fac_fill, ini_sol->c_done[c1]);
                FillFacility_2b(problem, ini_sol, fac_fill);
                cont++;
                break;
            }

        }

    }
    CompleteSolb(problem, ini_sol);
    for (int j = 0; j < problem->num_W; j++) {
        if (ini_sol->f_assigned[j] == 0 && ini_sol->fopen[j] == 1) {
            ini_sol->fopen[j] = 0;
            ini_sol->n_open--;
            ini_sol->cost = ini_sol->cost - problem->fixedcost[j];
        }
    }
}

void Heuristic::copy_solution(solution_t *sol1, const solution_t *sol2) {
    sol1->cost = sol2->cost;
    sol1->num_S = sol2->num_S;
    sol1->num_W = sol2->num_W;
    sol1->n_open = sol2->n_open;
    sol1->c_incompleted = sol2->c_incompleted;
    int S = sol2->num_S;
    int W = sol2->num_W;

    for (int i = 0; i < W; i++) {
        int longitud = sol2->facility_client[i][0];

            delete[] sol1->facility_client[i];
            sol1->facility_client[i] = new int[longitud + 1];

        std::memcpy(sol1->facility_client[i], sol2->facility_client[i], (longitud + 1) * sizeof(int));

        sol1->fopen[i] = sol2->fopen[i];
        sol1->f_assigned[i] = sol2->f_assigned[i];
    }

    for (int i = 0; i < S; i++) {
        int l = sol2->assignment_f[i][0];
            delete[] sol1->assignment_f[i];
            sol1->assignment_f[i] = new int[l + 1];
        std::memcpy(sol1->assignment_f[i], sol2->assignment_f[i], (l + 1) * sizeof(int));

            delete[] sol1->assignment_c[i];
            sol1->assignment_c[i] = new int[l + 1];
        std::memcpy(sol1->assignment_c[i], sol2->assignment_c[i], (l + 1) * sizeof(int));

        sol1->c_done[i] = sol2->c_done[i];
        sol1->g_assigned[i] = sol2->g_assigned[i];
    }

    for (int i = 0; i < S; i++) {
        int k = sol2->discards[i][1];

        delete[] sol1->discards[i];
        sol1->discards[i] = new int[k + 2];
        std::memcpy(sol1->discards[i], sol2->discards[i], (k + 2) * sizeof(int));

    }
}

void Heuristic::closef(const problem_t *problem, solution_t *close_sol, int facility) {
    int cliente, q, l2, pos2, estado_cli, f_open;
    int l = close_sol->facility_client[facility][0];
    f_open = close_sol->fopen[facility];
    vector<int> v_clientes;

        for (int i = 1; i < l+1; i++) {
            v_clientes.push_back(close_sol->facility_client[facility][i]);

        }

    while(!v_clientes.empty()) {

        cliente = v_clientes[0];
        estado_cli = close_sol->c_done[cliente];
        vector<int> info_cliente;
        vector<int> q_cliente;
        l2 = close_sol->assignment_f[cliente][0];


        for (int i = 0; i < l2+1; i++) {
            info_cliente.push_back(close_sol->assignment_f[cliente][i]);
            q_cliente.push_back(close_sol->assignment_c[cliente][i]);
        }

        auto pos = find(info_cliente.begin()+1, info_cliente.end(), facility);
        pos2 = distance(info_cliente.begin()+1,pos);
        q = q_cliente[pos2+1];
        close_sol->cost = close_sol->cost - q*problem->supplycost[cliente][facility];

        close_sol->f_assigned[facility] = close_sol->f_assigned[facility] - q;
        close_sol->g_assigned[cliente] = close_sol->g_assigned[cliente] - q;
        if (close_sol->g_assigned[cliente] < problem->goods[cliente]) {
            if (close_sol->g_assigned[cliente] == 0) {
                close_sol->c_done[cliente] = 0;
            }else {
                close_sol->c_done[cliente] = 1;
            }
        }
        if (estado_cli == 2 && close_sol->c_done[cliente] < 2) {
            close_sol->c_incompleted--;
        }

        info_cliente.erase(info_cliente.begin() + 1 + pos2);
        q_cliente.erase(q_cliente.begin() + 1 + pos2);

        info_cliente[0]--;
        q_cliente[0]--;
        delete [] close_sol->assignment_f[cliente];
        delete [] close_sol->assignment_c[cliente];

        close_sol->assignment_f[cliente] = new int[info_cliente[0]+1];
        close_sol->assignment_c[cliente] = new int[info_cliente[0]+1];
        for (int i = 0; i < info_cliente[0]+1; i++) {
            close_sol->assignment_f[cliente][i] = info_cliente[i];
            close_sol->assignment_c[cliente][i] = q_cliente[i];
        }

        v_clientes.erase(v_clientes.begin());
    }

    if (close_sol->f_assigned[facility] == 0) {
        close_sol->fopen[facility] = 0;
    }

    if (f_open == 1 && close_sol->fopen[facility] == 0) {
        close_sol->cost = close_sol->cost - problem->fixedcost[facility];

        close_sol->n_open--;

        delete [] close_sol->facility_client[facility];
        close_sol->facility_client[facility] = new int[1];

        close_sol->facility_client[facility][0] = 0;
    }

}

void Heuristic::ChangeSupplierb(const problem_t *problem, solution_t *solution, int cliente, int &millora) {

    if (solution->assignment_f[cliente][0] == 1) {

        int facility_old = solution->assignment_f[cliente][1];

        int pos_facility, facility_new;
        int c_goods = problem->goods[cliente];

        for (int i = 0; i < problem->num_W; i++) {
            if (problem->c_best_sc[cliente][i] == facility_old) {
                pos_facility = i;
                break;
            }
        }

        int coste_sc, coste_fixed;
        for (int j = 0; j < pos_facility; j++) {
            coste_sc = 0;
            coste_fixed = 0;
            facility_new = problem->c_best_sc[cliente][j];

            if (problem->supplycost[cliente][facility_new]<=Dmax) {
                int suma_inc = 0;
                int n_max = solution->facility_client[facility_new][0];
                int cliente2;


                for (int i = 1; i < n_max+1; i++) {
                    cliente2 = solution->facility_client[facility_new][i];
                    if (cliente != cliente2 && (problem->inc[cliente2][cliente] > 0 || problem->inc[cliente][cliente2] > 0)) {
                        suma_inc++;
                        break;
                    }

                }
                if (suma_inc == 0) {
                    if (problem->capacity[facility_new]-solution->f_assigned[facility_new]>=c_goods) {

                        coste_sc = c_goods*problem->supplycost[cliente][facility_new] - c_goods*problem->supplycost[cliente][facility_old];


                        if (solution->facility_client[facility_new][0]==0) {
                            coste_fixed = coste_fixed + problem->fixedcost[facility_new];
                        }
                        if (solution->facility_client[facility_old][0]==1) {
                            coste_fixed = coste_fixed - problem->fixedcost[facility_old];
                        }
                    }
                    if (coste_sc + coste_fixed < 0) {

                        DropCustomer(problem, solution, cliente, facility_old);
                        InsertCustomer_2(problem, solution, cliente, facility_new, solution->c_done[cliente]);
                        millora = 1;
                        break;
                    }
                }
            }

        }
    }else {
        vector<int> facility_pos;
        int max_pos;
        for (int i = 1; i < solution->assignment_f[cliente][0]+1; i++) {
            int facility_old = solution->assignment_f[cliente][i];

            for (int j = 0; j < problem->num_W; j++) {
                if (problem->c_best_sc[cliente][j] == facility_old) {
                    facility_pos.push_back(j);
                    if (i == 1) {
                        max_pos = j;
                    }else {
                        max_pos = maximo(max_pos, j);
                    }
                    break;
                }
            }
        }


        int facility_new;
        int coste_sc, coste_fixed;


        for (int j = 0; j < max_pos; j++) {
            coste_sc = 0;
            coste_fixed = 0;
            facility_new = problem->c_best_sc[cliente][j];

            if (problem->supplycost[cliente][facility_new]<=Dmax) {

                int suma_inc = 0;
                int n_max = solution->facility_client[facility_new][0];
                int cliente2;

                for (int i = 1; i < n_max+1; i++) {
                    cliente2 = solution->facility_client[facility_new][i];
                    if (cliente != cliente2 && (problem->inc[cliente2][cliente] > 0 || problem->inc[cliente][cliente2] > 0)) {
                        suma_inc++;
                        break;
                    }

                }

                int c_goods = problem->goods[cliente];
                if (suma_inc == 0) {
                    if (problem->capacity[facility_new]-solution->f_assigned[facility_new]>=c_goods) {

                        for (int k1 = 1; k1 < solution->assignment_f[cliente][0]+1; k1++) {
                            int fac = solution->assignment_f[cliente][k1];
                            coste_sc = coste_sc - solution->assignment_c[cliente][k1]*problem->supplycost[cliente][fac];

                            if (solution->facility_client[fac][0]==1) {
                                coste_fixed = coste_fixed - problem->fixedcost[fac];
                            }
                        }

                        coste_sc = coste_sc + c_goods*problem->supplycost[cliente][facility_new];

                        if (solution->facility_client[facility_new][0]==0) {

                            coste_fixed = coste_fixed + problem->fixedcost[facility_new];
                        }

                        if (coste_sc + coste_fixed < 0) {

                            for (int k = solution->assignment_f[cliente][0]; k > 0; k--) {
                                int facility_old = solution->assignment_f[cliente][k];
                                DropCustomer(problem, solution, cliente, facility_old);
                            }

                            InsertCustomer_2(problem, solution, cliente, facility_new, solution->c_done[cliente]);

                            millora = 1;
                            break;
                        }
                    }

                }
            }
        }
    }

}

void Heuristic::SplitCustomer3(const problem_t *problem, solution_t *solution, int cliente,int &millora) {

    millora = 0;
    int max_splits = max_cust;
    if (solution->assignment_f[cliente][0] < max_splits) {
        int max_pos = 0;
        int max_fac;
        int max_pos_assig = 0;
        vector<vector<int>> fac_pos;

        for (int i = 1; i < solution->assignment_f[cliente][0]+1; i++) {
            int facility_old = solution->assignment_f[cliente][i];

            for (int j = 0; j < problem->num_W; j++) {
                if (problem->c_best_sc[cliente][j] == facility_old) {
                    if (i == 1) {
                        max_pos = j;
                        max_fac = facility_old;
                        max_pos_assig = i;
                    }else {
                        if (j > max_pos) {
                            max_pos = maximo(j, max_pos);
                            max_fac = facility_old;
                            max_pos_assig = i;
                        }
                    }
                    fac_pos.push_back({facility_old,j,i});


                    break;
                }
            }
        }

        std::sort(fac_pos.begin(), fac_pos.end(), [](const std::vector<int>& a, const std::vector<int>& b) {
            return a[1] < b[1];
        });

        int fac_pos_size = fac_pos.size();

        int current_splits = solution->assignment_f[cliente][0];

        for(int i = 0; i < fac_pos_size; i++) {
            int facility = fac_pos[i][0];

            int position_fac = fac_pos[i][2];
            int q = solution->assignment_c[cliente][position_fac];
            if (cliente == 867) {

            }

            int max_pos_iter = fac_pos[i][1];
            int q1, q2, f1, f2;
            for (int f = 0; f < max_pos_iter; f++) {
                for (int g = f+1; g < max_pos_iter+1; g++) {
                    vector<vector<int>> new_matrix;

                    f1 = problem->c_best_sc[cliente][f];
                    f2 = problem->c_best_sc[cliente][g];

                    q1 = minimo(problem->capacity[f1]-solution->f_assigned[f1],q);
                    q2 = q - q1;

                    if (q2 <= problem->capacity[f2]-solution->f_assigned[f2]) {


                        if (q1 > 0 ) {
                            int suma_inc = 0;
                            int n_max = solution->facility_client[f1][0];
                            int cliente2;

                            if (problem->supplycost[cliente][f1]<=Dmax) {

                                for (int k = 1; k < n_max+1; k++) {

                                    cliente2 = solution->facility_client[f1][k];

                                    if (cliente != cliente2 && (problem->inc[cliente2][cliente] > 0 || problem->inc[cliente][cliente2] > 0)) {
                                        suma_inc++;
                                        break;
                                    }
                                }

                                if (suma_inc == 0) {


                                    new_matrix.push_back({f1,q1,q});
                                }
                            }


                        if (g != max_pos_iter && q2 != 0) {
                            if (problem->supplycost[cliente][f2]<=Dmax) {


                                suma_inc = 0;
                                n_max = solution->facility_client[f2][0];

                                for (int k = 1; k < n_max+1; k++) {

                                    cliente2 = solution->facility_client[f2][k];

                                    if (cliente != cliente2 && (problem->inc[cliente2][cliente] > 0 || problem->inc[cliente][cliente2] > 0)) {
                                        suma_inc++;
                                        break;
                                    }
                                }

                                if (suma_inc == 0 && q2 > 0) {
                                    new_matrix.push_back({f2,q2,q});
                                }
                            }
                        }

                        int m_size = new_matrix.size();

                        if (m_size > 0) {
                            int cost_sc = 0;
                            int cost_fixed = 0;

                            vector<int> fopen_temp(problem->num_W,0);

                            for (int rr = 0; rr < problem->num_W; rr++) {
                                fopen_temp[rr] = solution->fopen[rr];
                            }

                            for (int iter = 0; iter < m_size; iter++) {


                                int fac = new_matrix[iter][0];
                                int qtemp = new_matrix[iter][1];

                                cost_sc = cost_sc + qtemp*((int)problem->supplycost[cliente][fac]-(int)problem->supplycost[cliente][facility]);

                                if (fopen_temp[fac] == 0) {
                                    cost_fixed = cost_fixed + (int)problem->fixedcost[fac];
                                    fopen_temp[fac] = 1;
                                }
                            }

                            if ((m_size == 1 && new_matrix[0][1] == new_matrix[0][2] && solution->facility_client[facility][0] == m_size) || (m_size > 1 && new_matrix[0][2] == solution->f_assigned[facility])) {
                                cost_fixed = cost_fixed - (int)problem->fixedcost[facility];

                                fopen_temp[facility] = 0;
                            }

                            if (cost_sc + cost_fixed < 0) {

                                millora = 1;
                                current_splits++;
                                for (int iter = 0; iter < m_size; iter++) {
                                    if (cliente == 867) {
                                    }
                                    int fac2 = new_matrix[iter][0];
                                    int q3 = new_matrix[iter][1];
                                    int splits;

                                    if (m_size == 1) {
                                        if (new_matrix[iter][1]==new_matrix[iter][2]){
                                            splits = 1;
                                        }else {
                                            splits = 0;
                                        }
                                    }else {
                                        if (iter < m_size-1) {
                                            splits = 0;
                                        }else {
                                            splits = 1;
                                        }
                                    }

                                    DropCustomer_Split(problem, solution, cliente, facility, q3, splits);

                                    int done = solution->c_done[cliente];
                                    InsertCustomer_2q(problem, solution, cliente, fac2, q3, done);

                                }
                            }
                        }
                    }
                    }
                    if (millora == 1) {
                        break;
                    }
                }

                if (millora == 1) {
                    break;
                }
            }
            if (current_splits == max_splits || millora == 1) {
                break;
            }
        }


    }


}

vector<int> Heuristic::FacMenysCli(const problem_t *problem, solution_t *solution) {
    vector<vector<int>> numcli;

    for (int i = 0; i < problem->num_W; i++) {
        if (solution->fopen[i]==1) {
            numcli.push_back({i,solution->facility_client[i][0]});
        }
    }

    std::sort(numcli.begin(), numcli.end(), [](const std::vector<int>& a, const std::vector<int>& b) {
        return a[1] < b[1];
    });

    vector<int> FacNumCli;
    int max_size = numcli.size();
    for (int i = 0; i < max_size; i++) {
        FacNumCli.push_back(numcli[i][0]);

    }
    return FacNumCli;
}

vector<int> Heuristic::ClosestFac(const problem_t*problem, solution_t *solution, vector<int> v_cli, int facility_old) {

    int n = v_cli.size();
    int suma, cli, q;
    vector<vector<int>> SumCostFac;
    vector<int> vector_q;

    for (int i = 0; i < n; i++) {
        int cli = v_cli[i];
        int iter_max = solution->assignment_f[cli][0];

        for (int j = 1; j < iter_max+1; j++) {
            if (solution->assignment_f[cli][j] == facility_old) {
                vector_q.push_back(solution->assignment_c[cli][j]);
                break;
            }
        }

    }

    for (int j = 0; j < problem->num_W; j++) {

        if (solution->fopen[j]==0) {
            int cond = 0;
            for (int i = 0; i < n; i++) {
                q = vector_q[i];
                cli = v_cli[i];

                if (problem->supplycost[cli][j] <= Dmax) {
                    suma = suma + q*problem->supplycost[cli][j];
                }else {
                    cond = 1;
                    break;
                }

            }
            if (cond == 0) {
                suma = suma + problem->fixedcost[j];
                SumCostFac.push_back({j,suma});
            }
        }

    }

    std::sort(SumCostFac.begin(), SumCostFac.end(), [](const std::vector<int>& a, const std::vector<int>& b) {
        return a[1] < b[1];
    });


    vector<int> ReturnFac;

    for (int i = 0; i < SumCostFac.size(); i++) {
        ReturnFac.push_back(SumCostFac[i][0]);
    }

    return ReturnFac;
}

vector<int> Heuristic::FacMesCost(const problem_t *problem, solution_t *sol) {

    vector<vector<int>> RCL_fac;

    for (int i = 0; i < problem->num_W;i++) {
        if (sol->fopen[i]==1) {
            int n_cli = sol->facility_client[i][0];
            int suma = 0;

            for (int j = 1; j < n_cli+1; j++) {
                int cliente = sol->facility_client[i][j];
                int m = sol->assignment_f[cliente][0];
                for (int k = 1; k < m+1; k++) {
                    if (sol->assignment_f[cliente][k]== i) {
                        int q = sol->assignment_c[cliente][k];
                        suma = suma + q*problem->supplycost[cliente][i];
                        break;
                    }
                }

            }
            RCL_fac.push_back({i,suma});
        }
    }

    std::sort(RCL_fac.begin(), RCL_fac.end(), [](const std::vector<int>& a, const std::vector<int>& b) {
        return a[1] > b[1];
        });

    vector<int> resultado;
    for (const auto& fila : RCL_fac) {
        if (!fila.empty()) {
            resultado.push_back(fila[0]);
        }
    }
    return resultado;
}

void Heuristic::DesplazarInsertTabu(std::vector<int>& vec, int nuevoElemento, int modo) {
        if (modo == 0) {
            for (size_t i = vec.size() - 1; i > 0; --i) {
                vec[i] = vec[i - 1];
            }

            vec[0] = nuevoElemento;
        }else {
            auto it = std::find(vec.begin(), vec.end(), nuevoElemento);
            auto pos = distance(vec.begin(), it);

            for (size_t i = pos; i > 0; --i) {
                vec[i] = vec[i - 1];
            }
            vec[0] = nuevoElemento;
        }
}

vector<int> Heuristic::RatioFacCli(const problem_t * problem, solution_t * sol) {
    vector<vector<double>> RCL_fac;

    for (int i = 0; i < problem->num_W;i++) {
        if (sol->fopen[i]==1) {
            int n_cli = sol->facility_client[i][0];
            double suma = 0.0;

            for (int j = 1; j < n_cli+1; j++) {
                int cliente = sol->facility_client[i][j];
                int m = sol->assignment_f[cliente][0];
                for (int k = 1; k < m+1; k++) {
                    if (sol->assignment_f[cliente][k]== i) {
                        int q = sol->assignment_c[cliente][k];
                        suma = suma + q*problem->supplycost[cliente][i];
                        break;
                    }
                }

            }
            RCL_fac.push_back({(double)i,suma/n_cli});
        }
    }

    std::sort(RCL_fac.begin(), RCL_fac.end(), [](const std::vector<double>& a, const std::vector<double>& b) {
        return a[1] > b[1];
        });

    vector<int> resultado;
    for (const auto& fila : RCL_fac) {
        if (!fila.empty()) {
            resultado.push_back(fila[0]);
        }
    }
    return resultado;
}

void Heuristic::SwapIncompatibleCustomers(const problem_t* problem, solution_t *solution, int c1, int c2, int &millora) {
    int nfac1 = solution->assignment_f[c1][0];
    int nfac2 = solution->assignment_f[c2][0];

    int f1, f2, q1, q2;

    millora = 0;
    int exit_c = 0;

    for (int i = 1; i < nfac1+1; i++) {
        f1 = solution->assignment_f[c1][i];
        q1 = solution->assignment_c[c1][i];
        for (int j = 1; j < nfac2+1; j++) {

            f2 = solution->assignment_f[c2][j];
            q2 = solution->assignment_c[c2][j];

            if (problem->supplycost[c1][f2] <= Dmax) {

                if (problem->supplycost[c2][f1] <= Dmax) {
                    // Capacity f2 + q1 -q2
                    if (solution->f_assigned[f2] + q1 - q2 <= problem->capacity[f2]) {
                        // Capacity f1 + q2 - q1
                        if (solution->f_assigned[f1] + q2 - q1 <= problem->capacity[f1]) {
                            // incomp c1
                            int suma_inc = 0;
                            int n_max1 = solution->facility_client[f2][0];
                            int cliente2;


                            for (int k = 1; k < n_max1+1; k++) {
                                cliente2 = solution->facility_client[f2][k];

                                if (cliente2 != c2) {
                                    if (c1 != cliente2 && (problem->inc[cliente2][c1] > 0 || problem->inc[c1][cliente2] > 0)) {
                                        suma_inc++;
                                        break;
                                    }
                                }

                            }

                            if (suma_inc == 0) {
                                // incomp c2
                                int suma_inc2 = 0;
                                int n_max2 = solution->facility_client[f1][0];
                                int cliente2;


                                for (int p = 1; p < n_max2+1; p++) {
                                    cliente2 = solution->facility_client[f1][p];

                                    if (cliente2 != c1) {
                                        if (c2 != cliente2 && (problem->inc[cliente2][c2] > 0 || problem->inc[c2][cliente2] > 0)) {
                                            suma_inc2++;
                                            break;
                                        }
                                    }

                                }

                                if (suma_inc2 == 0) {

                                    int cost = problem->supplycost[c1][f2] + problem->supplycost[c2][f1] - problem->supplycost[c1][f1] - problem->supplycost[c2][f2];

                                    if (cost < 0) {
                                        millora = 1;

                                        // dropc1_f1
                                        DropCustomer(problem, solution, c1, f1);
                                        // dropc2_f2
                                        DropCustomer(problem, solution, c2, f2);
                                        // insertc1_f2
                                        InsertCustomer_2(problem, solution, c1, f2, solution->c_done[c1]);
                                        // insertc2_f1
                                        InsertCustomer_2(problem, solution, c2, f1, solution->c_done[c2]);

                                        exit_c = 1;
                                        break;

                                    }
                                }else{continue;}
                            }else{continue;}


                        }else{continue;}
                    }else{continue;}
                }else {continue;}
            }else {continue;}

        }

        if (exit_c == 1) {
            break;
        }
    }
}

void Heuristic::shakefinal(const problem_t *problem, solution_t *sol_t, solution_t *sol_bkp, solution_t *best, double &iter_shake, double kmins, double kmaxs,
    double ksteps, int &k, const std::chrono::steady_clock::time_point& start_time, int max_time, int &feasible) {

    if (sol_t->c_incompleted < problem->num_S) {
        cout << "break0.1"<<endl;
        exit(1);
    }
    copy_solution(sol_bkp,sol_t);

    int iter = 0;

    do {
        if (iter > 0) {
            copy_solution(sol_t,sol_bkp);
        }

        if (iter_shake == 1.0) {
            Constructive_8b(problem, sol_t);
            iter_shake = kmins;
        }
        else {
            int facopen = 0;
            for (int j = 0; j < problem->num_W; j++) {
                facopen = facopen + sol_t->fopen[j];
            }
            int num_pct_close = ceil(iter_shake*sol_t->n_open);
            int r = 0;

            while (r < num_pct_close) {
                int fac_close;
                int cond = 0;

                uniform_int_distribution<int> dist(0, problem->num_W-1);

                while (cond == 0) {
                    fac_close = dist(rng);
                    if (sol_t->fopen[fac_close]== 1) {
                        cond = 1;
                    }
                }

                closef(problem, sol_t, fac_close);

                r++;
            }

            CompleteSolb_shake(problem, sol_t);

        }
        if (sol_t->c_incompleted < problem->num_S) {
            iter++;
        }else {
            feasible = 1;
        }
        auto current_time = std::chrono::steady_clock::now();
        auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count();

        if (elapsed_time >= max_time) {
            break;
        }
    }while(feasible == 0);


    if (feasible == 1) {
        iter_shake = iter_shake + ksteps;
        k++;


        if (sol_t->cost < best->cost) {
            auto current_time = std::chrono::steady_clock::now();
            copy_solution(best,sol_t);
        }
    }

}


void Heuristic::vndfinal(const problem_t *problem, solution_t *sol_t, solution_t *best, int &k, const std::chrono::steady_clock::time_point& start_time, long long& best_time,
    long long elapsed_time, int max_time, int &break_cond, int &last_neigh) {

    while(k < max_iter2) {
                if (k == 4) {
                    if (sol_t->c_incompleted < problem->num_S) {
                        cout << "break2.1"<<endl;
                        exit(1);
                        break;
                    }

                    last_neigh = 4;
                    int cond;
                    vector<vector<int>> c_cost;
                    int cont_millores = 0;
                    vector<int> c_tabu(5,-1);

                    int millora;
                    do {
                        int maxpos;
                        cond = 0;

                        c_cost = CostPerCustomer(problem, sol_t);
                        vector<int> RCL_c_cost;

                        if (vnd_selec_cli == 2) {
                            double  pct = 0.1;
                            maxpos = ceil(c_cost.size()*pct);

                        }else if (vnd_selec_cli == 3){
                            int c_cost_size = c_cost.size();
                            int max = c_cost[0][1];
                            int min = c_cost[c_cost_size-1][1];
                            double limit_cost = max - 1.0*alpha*(max-min);

                            for (int r = 0; r < c_cost_size; r++) {
                                if (c_cost[r][1] >= limit_cost) {
                                    RCL_c_cost.push_back(c_cost[r][0]);
                                }else {
                                    break;
                                }
                            }
                            maxpos = RCL_c_cost.size();

                        }

                        uniform_int_distribution<int> dist_c_cost(0, maxpos-1);
                        millora = 0;

                        for (int j = 0; j < problem->num_S; j++) {

                            int cliente;

                            if (vnd_selec_cli == 1) {
                                cliente = c_cost[j][0];
                            }else if(vnd_selec_cli == 2) {
                                int temp_pos = dist_c_cost(rng);
                                cliente = c_cost[temp_pos][0];
                            }else if (vnd_selec_cli == 3){
                                int temp_pos = dist_c_cost(rng);
                                cliente = RCL_c_cost[temp_pos];
                            }

                            auto it = std::find(c_tabu.begin(), c_tabu.end(), cliente);

                            if (it == c_tabu.end()) {
                                ChangeSupplierb(problem, sol_t, cliente, millora);
                            }

                            if (millora == 1) {
                                DesplazarInsertTabu(c_tabu, cliente, 0);
                                cont_millores++;
                                break;
                            }

                            if (vnd_selec_cli == 2) {
                                if (j > 2*c_cost.size()) {
                                    cond = 1;
                                    break;
                                }
                            }else if (vnd_selec_cli == 3) {
                                if (j > 2*RCL_c_cost.size()) {
                                    cond = 1;
                                    break;
                                }
                            }
                        }

                        if (millora == 0) {
                            for (int j = 0; j < c_tabu.size(); j++) {
                                int cliente = c_tabu[j];
                                if (cliente != -1) {
                                    ChangeSupplierb(problem, sol_t, cliente, millora);
                                }

                                if (millora == 1) {
                                    DesplazarInsertTabu(c_tabu, cliente, 1);
                                    cont_millores++;
                                    cond = 0;
                                    break;
                                }
                                int max_j = c_tabu.size()-1;
                                if (j == max_j && millora == 0) {
                                    cond++;
                                }
                            }
                        }
                    }while(cond == 0);

                    if (k != 1) {
                        if (cont_millores == 0) {
                            k++;
                        }else {
                            k = 1;
                        }
                    }
                    else {
                        if (last_neigh == 2 && cont_millores == 0) {
                            k = k + 2;
                        }else {
                            k++;
                        }
                    }

                    if (sol_t->c_incompleted < problem->num_S) {
                        cout << "break2.2"<<endl;
                        exit(1);
                        break;
                    }
                    if (sol_t->cost < best->cost) {
                        auto current_time = std::chrono::steady_clock::now();
                        best_time = std::chrono::duration_cast<std::chrono::microseconds>(current_time - start_time).count();
                        copy_solution(best,sol_t);

                    }

                }
                else if (k == 3) {
                    if (sol_t->c_incompleted < problem->num_S) {
                        cout << "break3.1"<<endl;
                        exit(1);
                        break;
                    }

                    int c_ini = sol_t->cost;
                    last_neigh = 3;

                    int cond;
                    vector<vector<int>> c_cost;
                    int cont_millores = 0;

                    do {
                        cond = 0;
                        int maxpos;
                        c_cost = CostPerCustomer(problem, sol_t);

                        vector<int> RCL_c_cost;

                        if (vnd_selec_cli == 2) {
                            double  pct = 0.1;
                            maxpos = ceil(c_cost.size()*pct);
                        }else if (vnd_selec_cli == 3){
                            int c_cost_size = c_cost.size();
                            int max = c_cost[0][1];
                            int min = c_cost[c_cost_size-1][1];
                            double limit_cost = max - alpha*(max-min);


                            for (int r = 0; r < c_cost_size; r++) {
                                if (c_cost[r][1] >= limit_cost) {
                                    RCL_c_cost.push_back(c_cost[r][0]);
                                }else {
                                    break;
                                }
                            }
                            maxpos = RCL_c_cost.size();

                        }
                        uniform_int_distribution<int> dist_c_cost(0, maxpos-1);



                        int millora = 0;

                        for (int j = 0; j < problem->num_S; j++) {
                            int cliente;
                            if (vnd_selec_cli == 1) {
                                cliente = c_cost[j][0];
                            }else if(vnd_selec_cli == 2) {
                                int temp_pos = dist_c_cost(rng);
                                cliente = c_cost[temp_pos][0];
                            }else if (vnd_selec_cli == 3){
                                int temp_pos = dist_c_cost(rng);
                                cliente = RCL_c_cost[temp_pos];
                            }
                            int c1 = sol_t->cost;
                            SplitCustomer3(problem, sol_t,cliente,millora);
                           if (millora == 1) {
                                cont_millores++;

                            }
                            if (j == problem->num_S-1 && millora == 0) {
                                cond++;
                            }

                            if (vnd_selec_cli == 2) {
                                if (j > 2*c_cost.size()) {
                                    cond++;
                                    break;
                                }
                            }else if (vnd_selec_cli == 3) {
                                if (j > 2*RCL_c_cost.size()) {
                                    cond++;
                                    break;
                                }
                            }

                            auto current_time = std::chrono::steady_clock::now();
                            elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count();

                            if (millora == 1) {
                                break;
                            }else{
                                if (elapsed_time >= max_time) {
                                    break_cond = 1;
                                    break;
                                }
                            }

                        }

                        if (break_cond == 1) {
                            break;
                        }
                    }while(cond == 0);

                    if (k != 1) {
                        if (cont_millores == 0) {
                            k++;
                        }else {
                            k = 1;
                        }
                    }else {
                        if (last_neigh == 2 && cont_millores == 0) {
                            k = k + 2;
                        }else {
                            k++;
                        }
                    }

                    if (sol_t->c_incompleted < problem->num_S) {
                        cout << "break3.2"<<endl;
                        exit(1);
                        break;
                    }


                    if (sol_t->cost < best->cost) {
                        auto current_time = std::chrono::steady_clock::now();
                        best_time = std::chrono::duration_cast<std::chrono::microseconds>(current_time - start_time).count();
                        copy_solution(best,sol_t);

                    }
                    if (break_cond == 1) {
                        break;
                    }

                }
                else if (k == 2) {
                    if (sol_t->c_incompleted < problem->num_S) {
                        cout << "break4.1"<<endl;
                        exit(1);
                        break;
                    }
                    last_neigh = 2;
                    int cond = 0;
                    vector<vector<int>> c_cost;
                    int cont_millores = 0;
                    do {
                        int millora = 0;

                        vector<int> FacCli;

                        if (vnd_selec_fac == 1) {
                            FacCli = FacMenysCli(problem, sol_t);
                        }else if (vnd_selec_fac == 2) {
                            FacCli = FacMesCost(problem, sol_t);
                        }else {
                            FacCli = RatioFacCli(problem, sol_t);
                        }

                        int faccli_iter = FacCli.size();
                        int iter=0;

                        while(iter < faccli_iter) {
                            int maxpos;
                            int facility_close;
                            if (vnd_RCL_mode == 1) {
                                facility_close = FacCli[iter];
                            }else if (vnd_RCL_mode == 2) {
                                double  pct = 0.1;
                                maxpos = ceil(FacCli.size()*pct);

                            }
                            uniform_int_distribution<int> dist_selec_fac(0, maxpos-1);

                            if (vnd_RCL_mode == 2) {
                                int p = dist_selec_fac(rng);
                                facility_close = FacCli[p];
                            }

                            int ncli = sol_t->facility_client[facility_close][0];

                            vector<int> v_customers;
                            for (int j = 1; j < ncli+1; j++) {
                                v_customers.push_back(sol_t->facility_client[facility_close][j]);
                            }
                            int c_close = -(int)problem->fixedcost[facility_close];

                            vector<int> ClosestF = ClosestFac(problem, sol_t, v_customers, facility_close);
                            int facility_open = -1;

                            int n1 = ClosestF.size();
                            if (n1 > 0) {

                            for (int f = 0; f < n1; f++) {
                                int f_op = ClosestF[f];

                                if (f_op != facility_close && sol_t->fopen[f_op] == 0) {
                                        facility_open = f_op;
                                        break;
                                }
                            }
                            if (facility_open != -1) {
                                if (problem->capacity[facility_open] >= sol_t->f_assigned[facility_close]) {
                                    int c_open = (int)problem->fixedcost[facility_open];
                                    int coste_sc = 0;
                                    int coste_fixed = 0;

                                vector<vector<int>> new_matrix;
                                int estado, suma_estados;

                                vector<int> fopen_temp;

                                for (int f = 0; f < problem->num_W; f++) {
                                    fopen_temp.push_back(sol_t->fopen[f]);
                                }


                                int n2 = v_customers.size();
                                suma_estados = 0;


                                for(int r = 0; r < n2; r++) {

                                    int cliente = v_customers[r];

                                    int pos_f_close;

                                    for (int i2 = 1; i2 < sol_t->assignment_f[cliente][0]+1; i2++) {
                                        if (sol_t->assignment_f[cliente][i2] == facility_close) {
                                            pos_f_close = i2;
                                        }
                                    }

                                    int q = sol_t->assignment_c[cliente][pos_f_close];

                                     coste_sc = coste_sc + q*(problem->supplycost[cliente][facility_open]-problem->supplycost[cliente][facility_close]);

                                }

                                coste_fixed = c_open+c_close;
                                int c_initial = sol_t->cost;

                                if (coste_sc + coste_fixed < 0) {

                                    closef(problem, sol_t, facility_close);

                                    sol_t->cost = sol_t->cost + problem->fixedcost[facility_open];
                                    sol_t->fopen[facility_open] = 1;

                                    sol_t->n_open++;

                                    for (int c = 0; c < n2; c++) {
                                        int cli = v_customers[c];

                                        int done = sol_t->c_done[cli];
                                        InsertCustomer_2(problem, sol_t, cli, facility_open,done);
                                        }
                                    millora = 1;

                                    break;

                                }

                              }
                            }


                            }
                            iter++;
                        }


                        if (millora == 1) {
                            cont_millores++;
                        }
                        if (iter == faccli_iter && millora == 0) {
                            cond++;
                        }
                    }while(cond == 0);

                    int cafter = sol_t->cost;
                    if (k != 1) {
                        if (cont_millores == 0) {
                            k++;
                        }else {
                            k = 1;
                        }
                    }else {
                        if (last_neigh == 2 && cont_millores == 0) {
                            k = k + 2;
                        }else {
                            k++;
                        }
                    }
                    if (sol_t->c_incompleted < problem->num_S) {
                        cout << "break4.2"<<endl;
                        exit(1);
                        break;
                    }

                    if (sol_t->cost < best->cost) {
                        auto current_time = std::chrono::steady_clock::now();
                        best_time = std::chrono::duration_cast<std::chrono::microseconds>(current_time - start_time).count();
                        copy_solution(best,sol_t);
                    }
                }
                else if (k == 1){
                    if (sol_t->c_incompleted < problem->num_S) {
                        cout << "break5.1"<<endl;
                        exit(1);
                        break;
                    }

                    int cond = 0;
                    vector<vector<int>> c_cost;
                    int cont_millores = 0;

                    do {
                        int millora = 0;
                        vector<int> FacCli;
                        int maxpos;

                        if (vnd_selec_fac == 1) {
                            FacCli = FacMenysCli(problem, sol_t);
                        }else if (vnd_selec_fac == 2) {
                            FacCli = FacMesCost(problem, sol_t);
                            }else if (vnd_selec_fac == 3){
                            FacCli = RatioFacCli(problem, sol_t);
                            }

                        int faccli_iter = FacCli.size();


                        int iter=0;

                        while(iter < faccli_iter) {
                            int facility_close;
                            int maxpos;
                            if (vnd_RCL_mode == 1) {

                                facility_close = FacCli[iter];

                            }else if (vnd_RCL_mode == 2) {
                                double  pct = 0.1;
                                maxpos = ceil(FacCli.size()*pct);


                            }
                            uniform_int_distribution<int> dist_selec_fac(0, maxpos-1);

                            if (vnd_RCL_mode == 2) {
                                int p = dist_selec_fac(rng);
                                facility_close = FacCli[p];
                            }


                            int c_close = -(int)problem->fixedcost[facility_close];

                            int fac_pos;

                            int ncli = sol_t->facility_client[facility_close][0];

                            vector<vector<int>> v_customers_fpos;

                            for (int j = 1; j < ncli+1; j++) {
                                int c = sol_t->facility_client[facility_close][j];

                                int q;

                                for (int ii = 1; ii < sol_t->assignment_f[c][0]+1; ii++) {
                                    if (sol_t->assignment_f[c][ii] == facility_close) {
                                        q = sol_t->assignment_c[c][ii];
                                        break;
                                    }
                                }

                                for (int f = 0; f < problem->num_W; f++) {
                                    if (problem->c_best[c][f] == facility_close) {
                                        v_customers_fpos.push_back({c,f,q});
                                        break;
                                    }
                                }

                            }


                            int coste_sc = 0;
                            int estado;
                            int suma_estados = 0;
                            vector<int> fopen_temp;
                            vector<int> f_assigned_temp;
                            vector<vector<int>> new_matrix;
                            int m_size = 0;

                            for (int ff = 0; ff < problem->num_W; ff++) {
                                fopen_temp.push_back(sol_t->fopen[ff]);
                                f_assigned_temp.push_back(sol_t->f_assigned[ff]);

                            }

                            for (int i1 = 0; i1 < ncli; i1++) {
                                int c = v_customers_fpos[i1][0];


                                int max_pos = v_customers_fpos[i1][1];
                                int q = v_customers_fpos[i1][2];

                                if (max_pos != 0) {

                                    int fac_new = -1;
                                    for (int f1 = 0; f1 < max_pos; f1++) {
                                        fac_new = problem->c_best[c][f1];
                                        if (problem->supplycost[c][fac_new]<=Dmax) {
                                            if (fopen_temp[fac_new] == 1) {

                                                int suma_inc = 0;
                                                int n_max = sol_t->facility_client[fac_new][0];
                                                int cliente2;


                                                for (int k1 = 1; k1 < n_max+1; k1++) {
                                                    cliente2 = sol_t->facility_client[fac_new][k1];
                                                    if (c != cliente2 && (problem->inc[cliente2][c] > 0 || problem->inc[c][cliente2] > 0)) {
                                                        suma_inc++;
                                                        break;
                                                    }
                                                }

                                                if (suma_inc == 0 && problem->capacity[fac_new]-f_assigned_temp[fac_new]>=q) {


                                                    estado = 1;

                                                    coste_sc = coste_sc + q*(-problem->supplycost[c][facility_close] + problem->supplycost[c][fac_new]);


                                                    new_matrix.push_back({c,fac_new,estado});
                                                    f_assigned_temp[fac_new] = f_assigned_temp[fac_new] + q;
                                                    suma_estados = suma_estados + estado;
                                                    break;

                                                }else {
                                                    estado = 0;

                                                    new_matrix.push_back({c,-1,estado});
                                                    suma_estados = suma_estados + estado;
                                                    break;

                                                }
                                            }else if (f1 == max_pos -1){
                                                estado = 0;
                                                new_matrix.push_back({c,-1,estado});
                                            }
                                        }else if (f1 == max_pos -1){
                                            estado = 0;
                                            new_matrix.push_back({c,-1,estado});
                                            break;
                                        }
                                    }
                                    if (fac_new == -1) {
                                        new_matrix.push_back({c,-1,0});
                                        estado = 0;
                                        suma_estados = suma_estados + estado;

                                    }

                                }else {
                                    new_matrix.push_back({c,-1,0});
                                    estado = 0;
                                    suma_estados = suma_estados + estado;

                                }
                            }

                            m_size = new_matrix.size();

                            if (suma_estados == m_size && m_size > 0) {


                                int coste_total = coste_sc + c_close;


                                if (coste_total < 0) {

                                    closef(problem, sol_t, facility_close);
                                    for (int c1 = 0; c1 < new_matrix.size(); c1++) {

                                        int cli = new_matrix[c1][0];

                                        int done = sol_t->c_done[cli];
                                        InsertCustomer_2(problem, sol_t, new_matrix[c1][0], new_matrix[c1][1],done);
                                    }

                                    millora = 1;
                                    break;
                                }
                            }else {
                            }

                            iter++;
                        }

                        if (millora == 1) {
                            cont_millores++;
                        }
                        if (iter == faccli_iter && millora == 0) {
                            cond++;

                        }

                    }while(cond == 0);

                    if (k != 1) {
                        if (cont_millores == 0) {
                            k++;
                        }else {
                            k = 1;
                        }
                    }else {
                        if (last_neigh == 2 && cont_millores == 0) {
                            k = k + 2;
                        }else {
                            k++;
                        }
                    }

                    if (sol_t->c_incompleted < problem->num_S) {
                        cout << "break5.2"<<endl;
                        exit(1);
                    }

                    if (sol_t->cost < best->cost) {
                        auto current_time = std::chrono::steady_clock::now();
                        best_time = std::chrono::duration_cast<std::chrono::microseconds>(current_time - start_time).count();
                        copy_solution(best,sol_t);

                    }
                }
                else if (k == 5) {

                    if (sol_t->c_incompleted < problem->num_S) {
                        cout << "break5.1"<<endl;
                        exit(1);
                        break;
                    }
                    last_neigh = 5;

                    int cond;
                    vector<vector<int>> c_cost;
                    int cont_millores = 0;

                    vector<int> c_tabu(5,-1);
                    int millora;

                    do {
                        cond = 0;

                        c_cost = CostPerCustomerInc(problem, sol_t);

                         millora = 0;

                        for (int j = 0; j < c_cost.size(); j++) {

                            int cliente = c_cost[j][0];

                            vector<int> c_inc;

                            for (int r = 0; r < problem->num_S; r++) {
                                if (cliente != r) {
                                    if (problem->inc[cliente][r] == 1) {
                                        c_inc.push_back(r);
                                    }
                                }
                            }

                            if (c_inc.size() > 0) {
                                std::uniform_int_distribution<int> randomc(0, c_inc.size() - 1);

                                int randomIndex = randomc(rng);
                                int cliente2 = c_inc[randomIndex];

                                auto it = std::find(c_tabu.begin(), c_tabu.end(), cliente);

                                if (it == c_tabu.end()) {
                                    SwapIncompatibleCustomers(problem,sol_t,cliente,cliente2,millora);
                                }
                            }else{continue;}

                            if (millora == 1) {
                                DesplazarInsertTabu(c_tabu, cliente, 0);
                                cont_millores++;

                                break;
                            }

                        }

                        if (millora == 0) {
                            for (int j = 0; j < c_tabu.size(); j++) {

                                int cliente_tabu = c_tabu[j];

                                if (cliente_tabu != -1) {

                                    vector<int> c_inc;

                                    for (int r = 0; r < problem->num_S; r++) {
                                        if (cliente_tabu != r) {
                                            if (problem->inc[cliente_tabu][r] == 1) {
                                                c_inc.push_back(r);
                                            }
                                        }
                                    }

                                    if (c_inc.size() > 0) {
                                        std::uniform_int_distribution<int> randomc(0, c_inc.size() - 1);

                                        int randomIndex = randomc(rng);

                                        int cliente2 = c_inc[randomIndex];

                                        SwapIncompatibleCustomers(problem,sol_t,cliente_tabu,cliente2,millora);

                                    }else{continue;}

                                    if (millora == 1) {
                                        DesplazarInsertTabu(c_tabu, cliente_tabu, 0);
                                        cont_millores++;
                                        break;
                                    }


                                }

                                int max_j = c_tabu.size()-1;
                                if (j == max_j && millora == 0) {
                                    cond++;
                                }
                            }
                        }

                    }while(cond == 0);

                    if (k != 1) {
                        if (cont_millores == 0) {
                            k++;
                        }else {
                            k = 1;
                        }
                    }else {
                        if (last_neigh == 2 && cont_millores == 0) {
                            k = k + 2;
                        }else {
                            k++;
                        }
                    }

                    if (sol_t->c_incompleted < problem->num_S) {
                        cout << "break5.2"<<endl;
                        exit(1);
                        break;
                    }
                    if (sol_t->cost < best->cost) {
                        auto current_time = std::chrono::steady_clock::now();
                        best_time = std::chrono::duration_cast<std::chrono::microseconds>(current_time - start_time).count();
                        copy_solution(best,sol_t);

                    }

                }
            }
}

 void Heuristic::GVNSfinal(const problem_t *problem, solution_t *sol, long long elapsed_time, int max_time) {

    auto start_time = std::chrono::steady_clock::now();
    int k1;
    int i = 0;
    int cliente, facility, pos, pos_c;
    long long best_time = 0;
    uniform_int_distribution<int> dist(0, problem->num_S-1);
    int m1, m2, m3,m4, iter_best, millores_best, iter_best2, millores_best2;

    double const_bestsol;

    vector<int> tabu_fopen_ini(problem->num_W,0);
    auto best = create_solution(problem->num_S, problem->num_W);
    init_solution(problem, best);


    auto best_c = create_solution(problem->num_S, problem->num_W);
    init_solution(problem, best_c);

    auto sol_t = create_solution(problem->num_S, problem->num_W);
    init_solution(problem, sol_t);

    auto sol_bkp = create_solution(problem->num_S, problem->num_W);
    init_solution(problem, sol_bkp);

    long long elapsed_LS = 0;
    long long elapsed_Constr =0;
    millores_best2 = 0;
    iter_best2 = 0;
    int last_neigh;

    int n_feasible = 0;
    int n_infeasible = 0;

    auto start_Constr = std::chrono::steady_clock::now();

    int feas_cond = 0;

    do {
        Constructive_8b(problem, sol_t);
        if (sol_t->c_incompleted == problem->num_S) {
            feas_cond++;
        }
    }while(feas_cond == 0);

    auto current_Constr = std::chrono::steady_clock::now();
    elapsed_Constr = elapsed_Constr + std::chrono::duration_cast<std::chrono::microseconds>(current_Constr - start_Constr).count();
    double const_best = sol_t->cost;

    if (i == 0) {
        copy_solution(best,sol_t);
        copy_solution(best_c,sol_t);
        auto current_time = std::chrono::steady_clock::now();
        best_time = std::chrono::duration_cast<std::chrono::microseconds>(current_time - start_time).count();
        const_bestsol = const_best;
        i++;
    }

    int k = 1;
    int break_cond = 0;
    vndfinal(problem, sol_t, best, k, start_time, best_time, elapsed_time, max_time, break_cond, last_neigh);

    double ksteps = 1.0*(kmaxs-kmins)/nsteps;
    while (elapsed_time < max_time) {

        auto start_LS = std::chrono::steady_clock::now();

        millores_best = 0;
        iter_best = 0;

        double iter_shake = kmins;
        int feasible;
        while(iter_shake <= kmaxs) {
            int k = 0;
            int break_cond = 0;
            feasible = 0;
            shakefinal(problem, sol_t, sol_bkp, best, iter_shake, kmins, kmaxs, ksteps, k, start_time, max_time, feasible);

            if (feasible == 1) {
                vndfinal(problem, sol_t, best, k, start_time, best_time, elapsed_time, max_time, break_cond, last_neigh);
            }else {
                break;
            }


            if (break_cond == 1) {
                break;
            }
        }

        if(sol_t->c_incompleted < problem->num_S) {
            n_infeasible++;
            auto current_time = std::chrono::steady_clock::now();
            elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count();
            continue;
        }

        auto current_LS = std::chrono::steady_clock::now();
        elapsed_LS = elapsed_LS + std::chrono::duration_cast<std::chrono::microseconds>(current_LS - start_LS).count();
        auto current_time = std::chrono::steady_clock::now();
        elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count();

        n_feasible++;
    }


    copy_solution(sol, best);

    free_solution(best);
}
