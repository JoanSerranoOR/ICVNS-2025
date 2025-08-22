#ifndef GLOBALS_H
#define GLOBALS_H

#include <random>
using namespace std;

extern string name_inst;
extern int seed;
extern std::mt19937 rng;
extern double alpha;
extern int max_iter2;
extern double ini_open;
extern int c_best_mode;
extern int vnd_selec_cli;
extern int vnd_selec_fac;
extern int vnd_RCL_mode;
extern int Dmax;
extern int max_cust;
extern double kmins;
extern double kmaxs;
extern int nsteps;
#endif // GLOBALS_H