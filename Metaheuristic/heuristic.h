#ifndef HEURISTIC_H_
#define HEURISTIC_H_

#include "problem.h"
#include "solution.h"
#include <random>
#include <chrono>
using namespace std;

class Heuristic
{
public:
    Heuristic(problem_t *);
    ~Heuristic();

    problem_t *problem;

    void initialize();
    void init_solution(const problem_t *problem, solution_t *solution);
    void CompleteSolb(const problem_t *problem, solution_t *complete_sol);
    void CompleteSolb_shake(const problem_t *problem, solution_t *complete_sol);
    void Constructive_8b(const problem_t *problem, solution_t *ini_sol);
    vector<vector<int>> buildRCL1(const problem_t *problem, solution_t *sol);
    void InsertCustomer(const problem_t *problem, solution_t *sol, int cliente, int facility, int done);
    void DropCustomer(const problem_t *problem, solution_t *solution, int cliente, int facility_old);
    vector<vector<int>> CostPerCustomer(const problem_t *problem, const solution_t *sol);
    vector<vector<int>> CostPerCustomerInc(const problem_t *problem, const solution_t *sol);
    void InsertCustomer_2(const problem_t *problem, solution_t *sol, int cliente, int facility, int done);
    void InsertCustomer_2q(const problem_t *problem, solution_t *sol, int cliente, int facility, int suma_valor, int done);
    void FillFacility_2b(const problem_t *problem, solution_t *sol, int facility);
    void copy_solution(solution_t *sol1, const solution_t *sol2);
    void closef(const problem_t *problem, solution_t *solution, int facility);
    void ChangeSupplierb(const problem_t *problem, solution_t *solution, int cliente, int &millora);
    void SplitCustomer3(const problem_t *problem, solution_t *solution, int cliente,int &millora);
    void DropCustomer_Split(const problem_t *problem, solution_t *solution, int cliente, int facility_old, int q, int splits);
    vector<int> FacMenysCli(const problem_t *problem, solution_t *solution);
    vector<int> ClosestFac(const problem_t*problem, solution_t *solution, vector<int> v_cli, int facility_old);
    vector<int> RatioFacCli(const problem_t * problem, solution_t * sol);
    vector<int> FacMesCost(const problem_t *problem, solution_t *sol);
    void DesplazarInsertTabu(std::vector<int>& vec, int nuevoElemento, int modo);
    void GVNSfinal(const problem_t *problem, solution_t *sol, long long elapsed_time, int max_time);
    void SwapIncompatibleCustomers(const problem_t* problem, solution_t *solution, int c1, int c2, int &millora);
    void shakefinal(const problem_t *problem, solution_t *sol_t, solution_t *sol_bkp, solution_t *best, double &iter_shake, double kmins, double kmaxs,
    double ksteps, int &k, const std::chrono::steady_clock::time_point& start_time, int max_time, int &feasible);
    void vndfinal(const problem_t *problem, solution_t *sol_t, solution_t *best, int &k, const std::chrono::steady_clock::time_point& start_time, long long& best_time,
    long long elapsed_time, int max_time, int &break_cond, int &last_neigh);
};


#endif //HEURISTIC_H

