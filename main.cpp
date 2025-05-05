#include "Parameters.h"
#include "modelSetting.h"
#include "Lshape.h"
#include <chrono>
#include <set>

double prob_scenario(int s);
bool xi(int s, int k);
void solving(int dataNum);

int main() {
	cout << "Solve a subset of data files (press 0) or solve all data (press 1):";
	bool dataOption;
	cin >> dataOption;
	if (dataOption) {
		for (int dataNumber = 1; dataNumber <= 30; dataNumber++) {
			solving(dataNumber);
		}
	}
	else {
		cout << "How many data files?";
		int numFiles;
		cin >> numFiles;
		cout << "Input data numbers:";
		set<int> dataNumbers;
		int dataNumber;
		for (int i = 0; i < numFiles; i++) {
			cin >> dataNumber;
			dataNumbers.insert(dataNumber);
		}
		for (int dataNumber = 1; dataNumber <= 30; dataNumber++) {
			if(dataNumbers.find(dataNumber) != dataNumbers.end())
			solving(dataNumber);
		}
	}
	return 0;
}

void solving(int dataNum)
{
	//Elapsed time clock
	auto startTime = chrono::high_resolution_clock::now();

	//Basic setting: CPLEX setting and input parameters
	IloEnv env;
	inputParameters(env, "./data/data" + to_string(dataNum) + ".txt");
	cout << "\nStart solving data " << dataNum << "\n";

	//first-stage variables
	Var1D z(env, I_mobile, 0.0, 1.0), y(env, J_customer, 0.0, 1.0);
	Var2D w(env, I_mobile);
	for (int i = 0; i < I_mobile; i++) {
		w[i] = IloNumVarArray(env, J_customer, 0.0, 1.0);
	}
	Var1D theta = Var1D(env, sampleSize, 0.0, IloInfinity);

	//second-stage variables
	Var1D v(env, J_customer, 0.0, 1.0);
	Var3D p(env, I_mobile);
	for (int i = 0; i < I_mobile; i++) {
		p[i] = Var2D(env, J_customer);
		for (int j = 0; j < J_customer; j++) {
			p[i][j] = Var1D(env, K_shipper, 0.0, 1.0);
			for (int k = 0; k < K_shipper; k++) {
				//by constraint 9e
				if (t_mobile[i] != t_shipper[k]) p[i][j][k].setBounds(0.0, 0.0);
			}
		}
	}

	//set up first stage objective
	IloExpr objFirst = IloScalProd(cost_z, z) + IloScalProd(cost_y, y);
	double prob = 1.0 / sampleSize;
	for (int s = 0; s < sampleSize; s++) objFirst += prob * theta[s];

	//objective and constraints in second stage;
	IloExpr objSecond(env);
	IloRangeArray constr9b = IloRangeArray(env, K_shipper);
	IloRangeArray constr9c = IloRangeArray(env, J_customer);
	Constr2D constr9d = Constr2D(env, I_mobile);
	setSecondStage(env, objSecond, v, p, constr9b, constr9c, constr9d);

	//Optimal solution storage
	IloNum obj_opt;
	IloNumArray z_opt = IloNumArray(env, I_mobile);
	IloNumArray	y_opt = IloNumArray(env, J_customer);
	IloNumArray2D w_opt = IloNumArray2D(env, I_mobile);
	for (int i = 0; i < I_mobile; i++) {
		w_opt[i] = IloNumArray(env, J_customer);
	}
	IloNumArray theta_opt = IloNumArray(env, sampleSize);

	//Solution for estimation of upper bound
	IloNumArray z_upp(env, I_mobile);
	IloNumArray	y_upp(env, J_customer);
	IloNumArray2D w_upp(env, I_mobile);
	for (int i = 0; i < I_mobile; i++) {
		w_upp[i] = IloNumArray(env, J_customer);
	}

	/*Find the optimum*/
	IloNum realOptObj = 0;
	if (FindOptimum) {
		IloModel modelOpt(env);
		IloCplex cplexOpt(modelOpt);

		//(Optimum) first-stage variables
		IntVar1D z_real(env, I_mobile, 0.0, 1.0), y_real(env, J_customer, 0.0, 1.0);
		IntVar2D w_real(env, I_mobile);
		for (int i = 0; i < I_mobile; i++) {
			w_real[i] = IntVar1D(env, J_customer, 0.0, 1.0);
		}

		//(Optimum) second-stage variables
		int numScenario = pow(2, K_shipper);
		prob = 1.0 / numScenario;
		IntVar2D v_real(env, numScenario);
		IntVar4D p_real(env, numScenario);
		for (int s = 0; s < numScenario; s++) {
			v_real[s] = IntVar1D(env, J_customer, 0.0, 1.0);

			p_real[s] = IntVar3D(env, I_mobile);
			for (int i = 0; i < I_mobile; i++) {
				p_real[s][i] = IntVar2D(env, J_customer);
				for (int j = 0; j < J_customer; j++) {
					p_real[s][i][j] = IntVar1D(env, K_shipper, 0.0, 1.0);
					for (int k = 0; k < K_shipper; k++) {
						//by constraint 9e
						if (t_mobile[i] != t_shipper[k]) p_real[s][i][j][k].setBounds(0.0, 0.0);
					}
				}
			}
		}

		//(Optimum) objective
		IloExpr obj_real = IloScalProd(cost_z, z_real) + IloScalProd(cost_y, y_real);
		for (int s = 0; s < numScenario; s++) {
			for (int j = 0; j < J_customer; j++) {
				obj_real += (prob * cost_v[j]) * v_real[s][j];
			}
			for (int i = 0; i < I_mobile; i++) {
				for (int j = 0; j < J_customer; j++) {
					for (int k = 0; k < K_shipper; k++) {
						obj_real += (prob * cost_p[i][j][k]) * p_real[s][i][j][k];
					}
				}
			}
		}
		modelOpt.add(IloMinimize(env, obj_real));

		//constraits 7c
		for (int j = 0; j < J_customer; j++) {
			IloExpr sum(env);
			for (int i = 0; i < I_mobile; i++) {
				sum += w_real[i][j];
			}
			modelOpt.add(IloRange(sum + y_real[j] == 1.0));
			sum.end();
		}

		//constraints 7b and 7d
		for (int i = 0; i < I_mobile; i++) {
			modelOpt.add(IloRange(IloSum(w_real[i]) - C_capacity * z_real[i] <= 0.0));
		}

		
		for (int s = 0; s < numScenario; s++) {
			//constraints 7e
			for (int k = 0; k < K_shipper; k++) {
				IloExpr sum(env);
				for (int i = 0; i < I_mobile; i++) {
					for (int j = 0; j < J_customer; j++) {
						sum += p_real[s][i][j][k];
					}
				}
				modelOpt.add(IloRange(sum <= xi(s,k)));
				sum.end();
			}

			//constraint 7f
			for (int j = 0; j < J_customer; j++) {
				IloExpr sum(env);
				for (int i = 0; i < I_mobile; i++) sum += IloSum(p_real[s][i][j]);
				modelOpt.add(IloRange(sum + v_real[s][j] + y_real[j] == 1.0));
				sum.end();
			}

			//constraint 9d
			for (int i = 0; i < I_mobile; i++) {
				for (int j = 0; j < J_customer; j++) {
					modelOpt.add(IloRange(IloSum(p_real[s][i][j]) - w_real[i][j] <= 0.0));
				}
			}
		}

		cplexOpt.solve();

		ofstream output_opt("result_RealOptimum.txt");
		output_opt << cplexOpt.getObjValue() << "\n";
		IloNumArray z_result(env, I_mobile), y_result(env, J_customer);
		IloNumArray w_result(env, J_customer);
		cplexOpt.getValues(z_real, z_result);
		cplexOpt.getValues(y_real, y_result);
		output_opt << "z:" << z_result << endl;
		output_opt << "y:" << y_result << endl;
		for (int i = 0; i < I_mobile; i++) {
			cplexOpt.getValues(w_real[i], w_result);
			output_opt << "w" << i << ":" << w_result << endl;
		}
		output_opt.close();
	}


	//generate a sample to pick the sample for upper bound estimation
	IloArray<scenario> sample_pick(env, pickSize);
	for (int s = 0; s < pickSize; s++) sample_pick[s] = scenario(env, K_shipper);
	sampling(env, pickSize);
	assignSample(sample_pick);
	double min_U_mean = UINT32_MAX;

	vector<double> LowerBounds;
	vector<double> UpperBounds;

	for (int sampleCount = 0; sampleCount < numSamples; sampleCount++) {
		sample.end();
		sampling(env, sampleSize);
		LShape(env, objFirst, z, y, w, theta, v, p, obj_opt, z_opt, y_opt, w_opt, theta_opt, objSecond,
			constr9b, constr9c, constr9d);

		//change the right-hand side of constraints in second stage
		for (int j = 0; j < J_customer; j++) {
			constr9c[j].setBounds(1.0 - y_opt[j], 1.0 - y_opt[j]);
		}
		for (int i = 0; i < I_mobile; i++) {
			constr9d[i].setUbs(w_opt[i]);
		}

		//objective value (for lower bound estimation)
		obj_opt = IloScalProd(cost_z, z_opt) + IloScalProd(cost_y, y_opt);
		for (int s = 0; s < sampleSize; s++) {
			obj_opt += (1.0 / sampleSize) * secondStageObj(env, objSecond, constr9b, constr9c, constr9d, v, p, sample[s]);
		}
		//Lower bound
		double lowerbound = obj_opt;
		LowerBounds.push_back(lowerbound);

		//Upper bound
		obj_opt = IloScalProd(cost_z, z_opt) + IloScalProd(cost_y, y_opt);
		for (int s = 0; s < pickSize; s++) {
			obj_opt += (1.0 / pickSize) * secondStageObj(env, objSecond, constr9b, constr9c, constr9d, v, p, sample_pick[s]);
		}
		if (obj_opt < min_U_mean) {
			assignSol(z_opt, y_opt, w_opt, z_upp, y_upp, w_upp);
			min_U_mean = obj_opt;
		}
		cout << "Sample " << sampleCount + 1 << " finishes." << endl;
	}

	cout << "Start evaluating the upper bound.\n";
	//Evaluation of upper bound by sample_Eval
	//change the right-hand side of constraints in second stage
	for (int j = 0; j < J_customer; j++) {
		constr9c[j].setBounds(1.0 - y_upp[j], 1.0 - y_upp[j]);
	}
	for (int i = 0; i < I_mobile; i++) {
		constr9d[i].setUbs(w_upp[i]);
	}
	sample.end();
	sampling(env, evalSize);
	for (int s = 0; s < evalSize; s++) {
		double upper_bound = secondStageObj(env, objSecond, constr9b, constr9c, constr9d, v, p, sample[s]);
		UpperBounds.push_back(upper_bound);
		if (!(s % 1000) && s != 0) cout << s << " scenarios evaluated.\n";
	}

	//Statistics
	double mean_u = 0.0;
	double std_u = 0.0; //sample standard deviation
	double mean_l = 0.0;
	double std_l = 0.0; //sample standard deviation
	MeanAndStd(UpperBounds, mean_u, std_u);
	IloNum objValFirst_u = IloScalProd(cost_z, z_upp) + IloScalProd(cost_y, y_upp);
	mean_u += objValFirst_u;  //mean of upper bound must add the first stage objective value
	MeanAndStd(LowerBounds, mean_l, std_l);

	auto endTime = chrono::high_resolution_clock::now();
	auto elapsed = chrono::duration_cast<chrono::seconds>(endTime - startTime);

	ofstream output("./result/result" + to_string(dataNum) + ".txt");
	output << fixed << setprecision(4) << mean_u << " " << std_u << " ";
	output << mean_l << " " << std_l << " ";
	output << elapsed.count() << endl;
	output.close();

	cout << "Finished.\n";

	env.end();
}

double prob_scenario(int s)
{
	double p = 0;
	for (int k = 0; k < K_shipper; k++)
	{
		p *= s % 2 * p_cs + (1.0 - s % 2) * (1.0 - p_cs);
		s = s / 2;
	}
	return p;
}

bool xi(int s, int k)
{
	bool xi_value = 0;
	for (int d = 0; d < k; d++)
	{
		xi_value = s % 2;
		s = s / 2;
	}
	return xi_value;
}
