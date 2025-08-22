#include "globals.h"
#include <chrono>
using namespace std;

int seed = std::chrono::system_clock::now().time_since_epoch().count();
string name_inst = " ";

std::mt19937 rng(seed);
double alpha = 0.3; // 30
int max_iter2 = 6;
double ini_open = 0.3;
int c_best_mode = 3;//3
int vnd_selec_cli = 1;//1
int vnd_selec_fac = 1;
int vnd_RCL_mode = 1;//1
int Dmax = 0;
int max_cust = 2;
double kmins = 0.1;
double kmaxs = 1.0;
int nsteps = 1;