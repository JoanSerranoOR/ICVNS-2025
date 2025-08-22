#include "problem.h"
#include "heuristic.h"
#include "solution.h"
#include "globals.h"
//#include <QGridLayout>
//#include <QtCharts/QChartView>
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <fstream>
#include <random>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

problem_t *readinstance(problem_t *problem, string fichero);
void guardarSolucion3(solution_t *solucion, const string &file_name, int S, int W);
std::string getnameexec(const std::string& filename, const string num_iter);



// main irace -- define working directory previously
int main(int argc, char **argv) {
	problem_t *problem;
	string name_inst = argv[1];
	const string num_exec = argv[2];
	string name_sol_file = getnameexec(name_inst, num_exec);

	for (int i = 1; i < argc; ++i) {
		std::string arg_str = argv[i];

		if (arg_str == "--alpha") {
			if (i + 1 < argc) {
				try {
					alpha = std::stod(argv[++i]);
				} catch (const std::exception& e) {
					std::cerr << "ERROR: invalid value for --alpha: " << argv[i] << " (" << e.what() << ")" << std::endl;
					return 1;
				}
			} else {
				std::cerr << "ERROR: missing value for --alpha." << std::endl;
				return 1;
			}
		} else if (arg_str == "--ini_open") {
			if (i + 1 < argc) {
				try {
					ini_open = std::stod(argv[++i]);
				} catch (const std::exception& e) {
					std::cerr << "ERROR: invalid value for --ini_open: " << argv[i] << " (" << e.what() << ")" << std::endl;
					return 1;
				}
			} else {
				std::cerr << "ERROR: missing value for --ini_open." << std::endl;
				return 1;
			}
		}else if (arg_str == "--c_best_mode") {
			if (i + 1 < argc) {
				try {
					c_best_mode = std::atoi(argv[++i]);
				} catch (const std::exception& e) {
					std::cerr << "ERROR: invalid value for --c_best_mode: " << argv[i] << " (" << e.what() << ")" << std::endl;
					return 1;
				}
			} else {
				std::cerr << "ERROR: missing value for --c_best_mode." << std::endl;
				return 1;
			}
		}else if (arg_str == "--vnd_selec_fac") {
			if (i + 1 < argc) {
				try {
					vnd_selec_fac = std::atoi(argv[++i]);
				} catch (const std::exception& e) {
					std::cerr << "ERROR: invalid value for --vnd_selec_fac: " << argv[i] << " (" << e.what() << ")" << std::endl;
					return 1;
				}
			} else {
				std::cerr << "ERROR: missing value for --vnd_selec_fac." << std::endl;
				return 1;
			}
		}else if (arg_str == "--max_cust") {
			if (i + 1 < argc) {
				try {
					max_cust = std::atoi(argv[++i]);
				} catch (const std::exception& e) {
					std::cerr << "ERROR: invalid value for --max_cust: " << argv[i] << " (" << e.what() << ")" << std::endl;
					return 1;
				}
			} else {
				std::cerr << "ERROR: missing value for --max_cust." << std::endl;
				return 1;
			}
		}
		else if (arg_str == "--kmins") {
			if (i + 1 < argc) {
				try {
					kmins = std::stod(argv[++i]);
				} catch (const std::exception& e) {
					std::cerr << "ERROR: invalid value for --kmins: " << argv[i] << " (" << e.what() << ")" << std::endl;
					return 1;
				}
			} else {
				std::cerr << "ERROR: missing value for --kmins." << std::endl;
				return 1;
			}
		}else if (arg_str == "--kmaxs") {
			if (i + 1 < argc) {
				try {
					kmaxs = std::stod(argv[++i]);
				} catch (const std::exception& e) {
					std::cerr << "ERROR: invalid value for --kmaxs: " << argv[i] << " (" << e.what() << ")" << std::endl;
					return 1;
				}
			} else {
				std::cerr << "ERROR: missing value for --kmaxs." << std::endl;
				return 1;
			}
		}else if (arg_str == "--nsteps") {
			if (i + 1 < argc) {
				try {
					nsteps = std::atoi(argv[++i]);
				} catch (const std::exception& e) {
					std::cerr << "ERROR: invalid value for --nsteps: " << argv[i] << " (" << e.what() << ")" << std::endl;
					return 1;
				}
			} else {
				std::cerr << "ERROR: missing value for --nsteps." << std::endl;
				return 1;
			}
		}
	}

	problem = readinstance(problem, name_inst);
	const int time_limit1 = 10.0*sqrt(problem->num_W)*1000;

	Heuristic heur(problem);
	heur.initialize();

	auto *solution = create_solution(problem->num_S, problem->num_W);

	auto *best_solution = create_solution(problem->num_S, problem->num_W);
	heur.init_solution(problem,best_solution);

	long long elapsed_time1 = 0;

	heur.GVNSfinal(problem, solution,  elapsed_time1, time_limit1);
	cout << solution->cost << endl;
	guardarSolucion3(solution, "sol_" + name_sol_file, problem->num_S, problem->num_W);

	return 0;
}

problem_t *readinstance(problem_t *problem, const string fichero)
{
	int W, S;

	ifstream in(fichero);
	if (!in)
	{
		fprintf(stdout, "Error Instance Reading");
	}
	else {
		char *dummy = new char[200];
		char d;

		in.getline(dummy, 200, '=');
		in >> W;
		in.getline(dummy,200, '=');
		in >> S;


		problem = create_problem(S,W);

		in.getline(dummy, 200, '[');
		for (int j = 0; j < W; j++)
		{
			in >> problem->capacity[j] >> d;
		}


		in.getline(dummy, 200, '[');
		for (int j = 0; j < W; j++)
		{
			in >> problem->fixedcost[j] >> d;
		}

		in.getline(dummy, 200, '[');
		for (int j = 0; j < S; j++)
		{
			in >> problem->goods[j] >> d;

		}

		in.getline(dummy, 200, '[');
		in.getline(dummy, 200, '|');

		for (int i = 0; i < S; i++) {
			for (int j = 0; j < W; j++) {
				in >> problem->supplycost[i][j] >> d;
			}

		}



		in.getline(dummy,200,'=');
		in >> problem->num_Inc >> d;
		in.getline(dummy,200,'|');
		for (int i = 0; i < S; i++) {
			for (int ii = 0; ii < S; ii++) {
				problem->inc[i][ii] = 0;
			}
		}
		int fila,col;
		for (int i = 0; i < problem->num_Inc; i++) {
			in >> fila >> d;
			in >> col >> d;
			problem->inc[fila-1][col-1] = 1;
			problem->inc[col-1][fila-1] = 1;
		}

		in.getline(dummy,200,'=');
		in >> Dmax >> d;

	in.close();
    for (int i = 0; i < problem->num_W; i++) {
        auto **v_facilities = new int*[problem->num_S];
        for (int k = 0; k < problem->num_S; k++)
        {
            v_facilities[k] = new int[2];
        }
        for (int j = 0; j < problem->num_S; j++) {
            v_facilities[j][0] = j;
            v_facilities[j][1] = problem->supplycost[j][i];

        }

        std::sort(v_facilities, v_facilities + problem->num_S, [](const int* a, const int* b) {
            return a[1] < b[1];
        });
        for (int f = 0; f < problem->num_S; f++) {
            problem->f_best_c[i][f] = v_facilities[f][0];
        }
        for (int k = 0; k < problem->num_W; k++)
        {
            delete [] v_facilities[k];
        }
        delete [] v_facilities;
    }


	double min_scost;
	int facility;
	for (int i = 0; i < problem->num_S; i++) {
		if (c_best_mode == 1) {
			min_scost = problem->supplycost[i][0];

		}else if (c_best_mode == 2) {
			min_scost = problem->supplycost[i][0]+problem->fixedcost[0];
		}else {
			min_scost = problem->supplycost[i][0]+problem->fixedcost[0]/problem->capacity[0];
		}
		facility = 0;
		for (int j = 1; j < problem->num_W; j++) {
			if (c_best_mode == 1) {
				if (problem->supplycost[i][j] < min_scost) {
					min_scost = problem->supplycost[i][j];
					facility = j;
				}
			}else if (c_best_mode == 2) {
				if (problem->supplycost[i][j]+problem->fixedcost[j] < min_scost) {
					min_scost = problem->supplycost[i][j]+problem->fixedcost[j];
					facility = j;
				}
			}else {
				if (problem->supplycost[i][j]+problem->fixedcost[j]/problem->capacity[j] < min_scost) {
					min_scost = problem->supplycost[i][j]+problem->fixedcost[j]/problem->capacity[j];
					facility = j;
				}
			}


		}
		problem->c_best_f[i][0] = i;
		problem->c_best_f[i][1] = facility;
		problem->c_best_f[i][2] = min_scost;
	}

	std::sort(problem->c_best_f, problem->c_best_f + problem->num_S, [](const double* a, const double* b) {
		return a[2] > b[2];
	});


	for (int i = 0; i < problem->num_S; i++) {
		std::vector<std::vector<double>> v_facilities;
		std::vector<std::vector<double>> v_facilities_sc;

		for (int j = 0; j < problem->num_W; j++) {
			if (c_best_mode == 1) {
				v_facilities.push_back({(double)j,problem->supplycost[i][j]});
				v_facilities_sc.push_back({(double)j,problem->supplycost[i][j]});

			}else if (c_best_mode == 2) {
				v_facilities.push_back({(double)j,problem->supplycost[i][j]+problem->fixedcost[j]});
				v_facilities_sc.push_back({(double)j,problem->supplycost[i][j]});
			}else {
				v_facilities.push_back({(double)j,problem->supplycost[i][j]+problem->fixedcost[j]/problem->capacity[j]});
				v_facilities_sc.push_back({(double)j,problem->supplycost[i][j]});
			}

		}

		std::sort(v_facilities.begin(), v_facilities.end(), [](const std::vector<double>& a, const std::vector<double>& b) {
		return a[1] < b[1];
		});
		std::sort(v_facilities_sc.begin(), v_facilities_sc.end(), [](const std::vector<double>& a, const std::vector<double>& b) {
		return a[1] < b[1];
		});

		for (int k = 0; k < problem->num_W; k++) {
			problem->c_best[i][k] = (int)v_facilities[k][0];
			problem->c_best_sc[i][k] = (int)v_facilities_sc[k][0];
		}
	}

	int position;
	for (int i = 0; i < problem->num_S; i++) {
		problem->c_best_f_param[i][0] = i;

		for (int j = 0; j < problem->num_S; j++) {
			if (problem->c_best_f[j][0] == i) {
				position = j;
				break;
			}
		}
		problem->c_best_f_param[i][1] = position;
	}
	}
	return(problem);
}
void guardarSolucion3(solution_t *solucion, const string &file_name, int S, int W) {

	ofstream ofs(file_name);
	if (!ofs) {
		cerr << "Error opening " << file_name << "." << endl;
		return;
	}

	ofs << "[" << endl;
	int c;
	int copy;
	int copy2;
	for (unsigned i = 0; i < S; ++i) {
		ofs << "  (";
		for (unsigned j = 0; j < W; ++j) {
			c = 0;
			int l = solucion->assignment_f[i][0];
			for (unsigned k = 1; k < l+1; k++) {
				if (solucion->assignment_f[i][k] == j) {
					copy = solucion->assignment_c[i][k];
					c++;
				}else {
					copy2 = 0;
				}
			}
			if (c > 0) {
				ofs << copy;
			}else {
				ofs << copy2;
			}
			if (j < W - 1) {
				ofs << ", ";
			}
		}
		ofs << ")";
		ofs << endl;
	}
	ofs << "]" << endl;

	ofs.close();
}

std::string getnameexec(const std::string& filename, const string num_iter) {
	size_t lastDotTxtPos = filename.rfind(".txt");
	string return_nameinst;

	if (lastDotTxtPos != std::string::npos && lastDotTxtPos == filename.length() - 4) {
		return_nameinst = filename.substr(0, lastDotTxtPos) + "_" + num_iter + ".txt";
		return return_nameinst;
	} else {
		return filename;
	}
}