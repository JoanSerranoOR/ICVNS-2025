#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <gurobi_c++.h>
#include <fstream>
#include <cmath>
#include <map>
#include <regex>
#include <chrono>
#include <numeric>
using namespace std;

void read_data(const string fichero, vector<vector<int>> &supplycost, vector<int> &fixedcost, vector<int> &capacity,
	vector<int> &goods, vector<vector<int>> &incomp, int &W, int &S, int &ninc, int &Dmax) {
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

		cout << W << endl << S << endl;


		in.getline(dummy, 200, '[');
		int in_capacity;
		for (int j = 0; j < W; j++)
		{
			in >> in_capacity >> d;
			capacity.push_back(in_capacity);
		}


		in.getline(dummy, 200, '[');
		int in_fixed;
		for (int j = 0; j < W; j++)
		{
			in >> in_fixed >> d;
			fixedcost.push_back(in_fixed);
		}

		in.getline(dummy, 200, '[');
		int in_goods;
		for (int j = 0; j < S; j++)
		{
			in >> in_goods >> d;
			goods.push_back(in_goods);
		}


		in.getline(dummy, 200, '[');
		in.getline(dummy, 200, '|');


		supplycost = vector<vector<int>>(S, std::vector<int>(W, 0));
		int in_supp;
		for (int i = 0; i < S; i++) {
			for (int j = 0; j < W; j++) {
				in >> in_supp >> d;
				supplycost[i][j] = in_supp;
			}

		}

		in.getline(dummy,200,'=');
		in >> ninc >> d;
		in.getline(dummy,200,'|');

		incomp = vector<vector<int>>(S, std::vector<int>(S, 0));

		int fila,col;
		for (int i = 0; i < ninc; i++) {
			in >> fila >> d;
			in >> col >> d;
			incomp[fila-1][col-1] = 1;
			incomp[col-1][fila-1] = 1;
		}

		in.getline(dummy, 200, '=');
		in >> Dmax;

		in.close();
	}
}

void addConstraintsToGurobiModel2(GRBModel& model,
                                 std::map<std::pair<int,int>, GRBVar> X,
                                 GRBVar* Y,
                                 std::map<std::pair<int,int>, GRBVar> Z,
                                 vector<int> goods,
                                 vector<int> capacities,
                                 vector<vector<int>> incomp,
                                 vector<vector<int>> supplycost,
                                 int S_val, int W_val, int Dmax)
{
    try {
        for (int j = 0; j < W_val; ++j) { 
            GRBLinExpr lhs = 0; 

              for (int i = 0; i < S_val; ++i) { 
            	if (X.count({i, j})) {
            		lhs += X[{i,j}]; 
            	}
            }

            
            lhs -= capacities[j] * Y[j]; 

            std::string constraint_name = "Capacity_" + std::to_string(j);
            model.addConstr(lhs <= 0, constraint_name);
            
        }

        for (int i = 0; i < S_val; ++i) {
            GRBLinExpr lhs = 0;

            for (int j = 0; j < W_val; ++j) {
            	if (X.count({i, j})) {
            		lhs += X[{i,j}];
            	}
            }

            std::string constraint_name = "Goods_" + std::to_string(i);
            model.addConstr(lhs == goods[i], constraint_name);
        }

    	GRBLinExpr lhs;
    	for (int i = 0; i < S_val; ++i) { 
    		for (int j = 0; j < W_val; ++j) { 

    			if (X.count({i, j})) {
    				lhs = X[{i,j}] - goods[i]*Z[{i,j}];
    			}

    			std::string constraint_name = "UpperX_" + std::to_string(i) + "_" + std::to_string(j);
    			model.addConstr(lhs <= 0, constraint_name);
    		}

    	}

    	  for (int i1 = 0; i1 < S_val; ++i1) { 
            for (int i2 = i1 + 1; i2 < S_val; ++i2) { 

                if (incomp[i1][i2] == 1) {

                    for (int j = 0; j < W_val; ++j) {
                    	if (Z.count({i1, j}) && Z.count({i2, j})) {
                    		GRBLinExpr lhs = Z[{i1,j}] + Z[{i2,j}];
                    		std::string constr_name1 = "Incomp_" + std::to_string(i1) + "_" + std::to_string(i2) + "_" + std::to_string(j);
                    		model.addConstr(lhs <= Y[j], constr_name1);
                    	}
                    }



                }
            }
        }




        std::cout << "Restricciones añadidas exitosamente." << std::endl;

    } catch (GRBException e) {
        std::cerr << "Error Gurobi: " << e.getMessage() << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

void addConstraintsToGurobiModel3(GRBModel& model,
                                 GRBVar** X,
                                 GRBVar** Z,
                                 GRBVar* Y,
                                 vector<int> goods,
                                 vector<int> capacities,
                                 vector<vector<int>> incomp,
                                 vector<vector<int>> supplycost,
                                 int S_val, int W_val)
{
    try {
        for (int j = 0; j < W_val; ++j) {
            GRBLinExpr lhs = 0;

            for (int i = 0; i < S_val; ++i) {
            		lhs += X[i][j];
            }

            lhs -= capacities[j] * Y[j];

            std::string constraint_name = "Capacity_" + std::to_string(j);
            model.addConstr(lhs <= 0, constraint_name);
        }

        for (int i = 0; i < S_val; ++i) {
            GRBLinExpr lhs = 0;

            for (int j = 0; j < W_val; ++j) {
            		lhs += X[i][j];
            }

            std::string constraint_name = "Goods_" + std::to_string(i);
            model.addConstr(lhs == goods[i], constraint_name);
        }

    	GRBLinExpr lhs;
    	for (int i = 0; i < S_val; ++i) {
    		for (int j = 0; j < W_val; ++j) {
    			lhs = X[i][j] - goods[i]*Z[i][j];


    			std::string constraint_name = "UpperX_" + std::to_string(i) + "_" + std::to_string(j);
    			model.addConstr(lhs <= 0, constraint_name);
    		}

    	}


    	for (int i1 = 0; i1 < S_val; ++i1) {
    		for (int i2 = i1 + 1; i2 < S_val; ++i2) {
    			if (incomp[i1][i2] == 1) {

    				for (int j = 0; j < W_val; ++j) {
    						GRBLinExpr lhs = Z[i1][j] + Z[i2][j];
    						model.addConstr(lhs <= Y[j]); 
    				}


    			}
    		}
    	}

        std::cout << "Constraints Added." << std::endl;

    } catch (GRBException e) {
        std::cerr << "Error Gurobi: " << e.getMessage() << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

void displayGurobiSolution(GRBModel& model, std::map<std::pair<int,int>, GRBVar> X, GRBVar* Y, int S, int W) {
    try {
        int status = model.get(GRB_IntAttr_Status);

        if (status == GRB_OPTIMAL) {
            std::cout << "Solucion optima encontrada:" << std::endl;

            bool all_integers = true;
            double tolerance = 1e-6;

            for (int i = 0; i < S; ++i) {
                for (int j = 0; j < W; ++j) {

                	if (X.count({i, j})) {
                		double val_X = X[{i,j}].get(GRB_DoubleAttr_X);

                		if (std::abs(val_X - std::round(val_X)) > tolerance) {
                			all_integers = false;
                			break;
                		}
                	}


                }
                if (!all_integers) break;
            }


            if (all_integers) {
                for (int j = 0; j < W; ++j) {
                    double val_Y = Y[j].get(GRB_DoubleAttr_X);
                    if (std::abs(val_Y - 0.0) > tolerance && std::abs(val_Y - 1.0) > tolerance) {
                        all_integers = false;
                        break;
                    }
                }
            }

            if (all_integers) {
                std::cout << "Integer Solution Found." << std::endl;
            } else {
                 std::cout << "Optimal not integer solution found." << std::endl;
            }

            std::cout << "\nObj Value = " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;



        } 
    } catch (GRBException e) {
        std::cerr << "Error Code: " << e.getErrorCode() << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Standard Error: " << e.what() << std::endl;
    }
}

void displayGurobiSolution3(GRBModel& model, GRBVar** X, GRBVar* Y, int S, int W) {
    try {
        int status = model.get(GRB_IntAttr_Status);

        if (status == GRB_OPTIMAL) {
            std::cout << "Solucion optima encontrada:" << std::endl;

            bool all_integers = true;
            double tolerance = 1e-6;

            for (int i = 0; i < S; ++i) {
                for (int j = 0; j < W; ++j) {

                		double val_X = X[i][j].get(GRB_DoubleAttr_X);
                		if (std::abs(val_X - std::round(val_X)) > tolerance) {
                			all_integers = false;
                			break;
                		}

                }
                if (!all_integers) break;
            }

            if (all_integers) {
                for (int j = 0; j < W; ++j) {
                    double val_Y = Y[j].get(GRB_DoubleAttr_X);

                    if (std::abs(val_Y - 0.0) > tolerance && std::abs(val_Y - 1.0) > tolerance) {
                        all_integers = false;
                        break;
                    }
                }
            }

            if (all_integers) {
                std::cout << "Integer Solution Found." << std::endl;
            } else {
                std::cout << "Optimal not integer Solution Found,." << std::endl;
            }

            std::cout << "\nObj Value = " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;



        } 
    }catch (GRBException e) {
        std::cerr << "Error code: " << e.getErrorCode() << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Standard Error : " << e.what() << std::endl;
    }
}

void saveToFile(const std::string &filename, std::vector<double> data) {
	std::ofstream outFile(filename, std::ios::app);
	if (outFile.is_open()) {
		outFile << std::fixed; 

		for (const auto &value : data) {
			std::ostringstream oss;
			oss << value;

			std::string strValue = oss.str();
			std::replace(strValue.begin(), strValue.end(), '.', ',');

			outFile << strValue << " ";
		}

		outFile << "\n"; 
		outFile.close();
	} else {
		std::cerr << "The file could not be opened." << std::endl;
	}
}


std::string getnameinstMD(const std::string& filename) {
	size_t lastDotTxtPos = filename.rfind(".txt");
	string return_nameinst;

	if (lastDotTxtPos != std::string::npos && lastDotTxtPos == filename.length() - 4) {
		return_nameinst = filename.substr(0, lastDotTxtPos) + "_MB.txt";
		return return_nameinst;
	} else {
		return filename;
	}
}

int main(int argc, char **argv) {
    try {

        vector<vector<int>> supplycost;
        vector<int> fixedcost;
        vector<int> capacity;
        vector<int> goods;
        vector<vector<int>> incomp;
    	int n_sup, n_cust, ninc, Dmax, DmaxIni;

		string name_inst = argv[1];
    	int exact_mode = atoi(argv[2]);
    	

    	string name_instMD = getnameinstMD(name_inst);
    	read_data(name_instMD, supplycost, fixedcost, capacity, goods, incomp, n_sup, n_cust, ninc, DmaxIni);


    	vector<vector<int>> ini_x = vector(n_cust, std::vector<int>(n_sup, 0));
    	vector<int> ini_y = vector(n_sup, 0);

    
        GRBEnv env = GRBEnv();
        GRBModel model = GRBModel(env);

        model.set(GRB_StringAttr_ModelName, "MyModel");

    	GRBVar** X1;
    	GRBVar** Z1;
    	std::map<std::pair<int,int>, GRBVar> X2;
    	std::map<std::pair<int,int>, GRBVar> Z2;
    	if (exact_mode == 1 || exact_mode == 3) {

    		X1 = new GRBVar*[n_cust];
    		Z1 = new GRBVar*[n_cust];

    		for (int i = 0; i < n_cust; ++i) {
    			X1[i] = new GRBVar[n_sup];
    			Z1[i] = new GRBVar[n_sup];
    			for (int j = 0; j < n_sup; ++j) {

    				std::string var_name = "X_" + std::to_string(i) + "_" + std::to_string(j);
    				X1[i][j] = model.addVar(0.0, goods[i], 0.0, GRB_INTEGER, var_name);
    				std::string var_name2 = "Z_" + std::to_string(i) + "_" + std::to_string(j);
    				Z1[i][j] = model.addVar(0.0, 1, 0.0, GRB_BINARY, var_name2);
    			}
    		}
    	}else {
    		for (int i = 0; i < n_cust; ++i) {
    			for (int j = 0; j < n_sup; ++j) {

    				if (supplycost[i][j] > DmaxIni) {
    					continue;
    				}

    				std::string var_name = "X_" + std::to_string(i) + "_" + std::to_string(j);
    				X2[{i, j}] = model.addVar(
						0.0, goods[i], 0.0, GRB_INTEGER, var_name
					);

    				std::string var_name2 = "Z_" + std::to_string(i) + "_" + std::to_string(j);
    				Z2[{i,j}] = model.addVar(0.0, 1, 0.0, GRB_BINARY, var_name2);

    			}
    		}
    	}

    	GRBVar* Y = new GRBVar[n_sup];
    	for (int j = 0; j < n_sup; ++j) {
    		std::string var_name = "Y_" + std::to_string(j);
    		Y[j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
    		Y[j].set(GRB_DoubleAttr_Start, ini_y[j]);
    	}

    	GRBLinExpr obj = 0;

    	if (exact_mode == 1 || exact_mode == 3) {
    		for (int i = 0; i < n_cust; ++i) {
    			for (int j = 0; j < n_sup; ++j) {
    				obj += static_cast<double>(supplycost[i][j]) * X1[i][j];
    			}
    		}
    	}else {
    		for (int i = 0; i < n_cust; ++i) {
    			for (int j = 0; j < n_sup; ++j) {

    				if (X2.count({i, j})) {
    					obj += static_cast<double>(supplycost[i][j]) * X2[{i,j}];
    				}

    			}
    		}
    	}


    	for (int j = 0; j < n_sup; ++j) {
    		obj += static_cast<double>(fixedcost[j]) * Y[j];
    	}

    	model.setObjective(obj, GRB_MINIMIZE);


    	cout << "pre restricciones"<< endl;
        if (exact_mode == 1 || exact_mode == 3) {
    		addConstraintsToGurobiModel3(model,X1,Z1,Y,goods,capacity, incomp, supplycost,n_cust,n_sup);
    	}else {
    		addConstraintsToGurobiModel2(model,X2,Y,Z2,goods,capacity, incomp, supplycost,n_cust,n_sup, DmaxIni);
    	}


		if (exact_mode == 3) {
			// ini sol
			std::vector<int> valores_iniciales_y_vec(n_sup, 0);
			std::vector<std::vector<int>> valores_iniciales_x_mat(n_cust, std::vector<int>(n_sup, 0));

		

			for (int i = 0; i < n_cust; ++i) {
				for (int j = 0; j < n_sup; ++j) {
					X1[i][j].set(GRB_DoubleAttr_Start, valores_iniciales_x_mat[i][j]);

					if (valores_iniciales_x_mat[i][j]>0) {
						Z1[i][j].set(GRB_DoubleAttr_Start, 1);
					}else {
						Z1[i][j].set(GRB_DoubleAttr_Start, 0);
					}
					}
			}


			for (int j = 0; j < n_sup; ++j) {
				for (int i = 0; i < n_cust; ++i) {
					if (valores_iniciales_x_mat[i][j]>0) {
						Y[j].set(GRB_DoubleAttr_Start, 1);
						break;
					}

					}
			}

		}

    	model.update();

    	model.set(GRB_DoubleParam_TimeLimit, 3600);
    	model.set(GRB_IntParam_OutputFlag, 1);
    	model.set(GRB_DoubleParam_MIPGap, 0.0);

    	std::string baseFileName;

    	size_t lastDotPos = name_inst.rfind('.');

    	if (lastDotPos != std::string::npos && lastDotPos < name_inst.length() - 1) {
    		baseFileName = name_inst.substr(0, lastDotPos);
    	} else {
    		baseFileName = name_inst;
    	}

    	model.set(GRB_StringParam_LogFile, "gurobi_output_" + baseFileName + "_" + to_string(DmaxIni) + ".log"); 
    	model.set(GRB_IntParam_LogToConsole, 0);

    	auto start_time = std::chrono::high_resolution_clock::now();

        model.optimize();

    	auto end_time = std::chrono::high_resolution_clock::now();
    	std::chrono::duration<double> duration = end_time - start_time;

    	double runtime = model.get(GRB_DoubleAttr_Runtime);

    	if (exact_mode == 1) {
    		displayGurobiSolution3(model, X1, Y, n_cust, n_sup);
    	}else {
    		displayGurobiSolution(model, X2, Y, n_cust, n_sup);
    	}

    	/* to save the results
    	 *std::vector<double> data;
    	if (model.get(GRB_IntAttr_SolCount) > 0) {
    		data = {model.get(GRB_DoubleAttr_ObjVal),model.get(GRB_DoubleAttr_ObjBound),model.get(GRB_DoubleAttr_MIPGap),runtime, 0.0,(double)DmaxIni,	(double)Dmax};
    	}else {
    		data = {0.0,0.0,0.0, (double)DmaxIni,	(double)Dmax};
    	}

		saveToFile("Results.txt",data);*/

    } catch (GRBException e) {
        std::cerr << "Error Gurobi: " << e.getMessage() << std::endl;
        std::cerr << "Error Code: " << e.getErrorCode() << std::endl;
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }


    return 0;
}
