#include "Parameters.h"
#include <fstream>
#include <unordered_set>
#include <numeric>

//Sets and parameters
bool FindOptimum;
set_cardinality I_mobile, J_customer, K_shipper;
IloNumArray cost_z;
IloNumArray3D cost_p;
IloNumArray cost_y;
IloNumArray cost_v;
IloBoolArray t_mobile;
IloBoolArray t_shipper;
int sampleSize;
int numSamples;
int pickSize;
int evalSize;
double p_cs;
int C_capacity;
IloArray<scenario> sample;

double** Allocate2DArray(int len1, int len2);
void delete2DArray(double** array, int len1);

void inputParameters(IloEnv env, string filename) {
	
	ifstream fin(filename);

	//Sample parameters
	fin >> FindOptimum;
	fin >> sampleSize;
	fin >> numSamples;
	fin >> pickSize;
	fin >> evalSize;

	//Cardinalities
	fin >> I_mobile;
	fin >> J_customer;
	fin >> K_shipper;
	
	double** mobile = Allocate2DArray(I_mobile, 2);
	double** customer = Allocate2DArray(J_customer, 2);
	double** shipper = Allocate2DArray(K_shipper, 4);

	cost_z = IloNumArray(env, I_mobile);
	cost_p = IloNumArray3D(env, I_mobile);
	for (int i = 0; i < I_mobile; i++) {
		cost_p[i] = IloNumArray2D(env, J_customer);
		for(int j=0;j<J_customer;j++) cost_p[i][j] = IloNumArray(env, K_shipper);
	}

	cost_y = IloNumArray(env, J_customer);
	cost_v = IloNumArray(env, J_customer);
	t_mobile = IloBoolArray(env, I_mobile);
	t_shipper = IloBoolArray(env, K_shipper);

	//Position of mobile depots
	cout << "\nPosition of mobile depots:\n";
	for (int i = 0; i < I_mobile; i++)
	{
		cout << i << "\t";
		for (int x = 0; x < 2; x++)
		{
			fin >> mobile[i][x];
			cout << mobile[i][x] << "\t";
		}
		cout << "\n";
	}
	cout << "\n";

	//Position of customers
	cout << "\nPosition of customers:\n";
	for (int j = 0; j < J_customer; j++)
	{
		cout << j << "\t";
		for (int x = 0; x < 2; x++)
		{
			fin >> customer[j][x];
			cout << customer[j][x] << "\t";
		}
		cout << "\n";
	}
	cout << "\n";

	//Position of crowd-shippers
	cout << "\nPosition of crowd-shippers:\n";
	for (int k = 0; k < K_shipper; k++)
	{
		cout << k << "\t";
		for (int x = 0; x < 4; x++)
		{
			fin >> shipper[k][x];
			cout << shipper[k][x] << "\t";
		}
		cout << "\n";
	}

	//Time availability of mobile depots
	cout << "\nTime available (mobile depots):\n";
	for (int i = 0; i < I_mobile; i++)
	{
		fin >> t_mobile[i];
		cout << t_mobile[i] << "\t";
	}

	//Time availability of crowd-shippers
	cout << "\nTime available (crowd-shippers):\n";
	for (int k = 0; k < K_shipper; k++)
	{
		fin >> t_shipper[k];
		cout << t_shipper[k] << "\t";
	}
	cout << "\n";

	//Probability of crowd-shipper availability
	fin >> p_cs;
	//Capacity
	fin >> C_capacity;
	fin.close();

	//Costs
	cout << "cost_z\n";
	for (int i = 0; i < I_mobile; i++)
	{
		cost_z[i] = sqrt(pow(mobile[i][0], 2) + pow(mobile[i][1], 2));
		cout << cost_z[i] << "\t";
	}
	cout << "\n";

	cout << "cost_p\n";
	for (int i = 0; i < I_mobile; i++)
	{
		for (int j = 0; j < J_customer; j++)
		{
			for (int k = 0; k < K_shipper; k++)
			{
				cost_p[i][j][k] = 1.0 + 0.6 * \
					(sqrt(pow(mobile[i][0] - shipper[k][0], 2) + pow(mobile[i][1] - shipper[k][1], 2))\
						+ sqrt(pow(customer[j][0] - mobile[i][0], 2) + pow(customer[j][1] - mobile[i][1], 2))\
						+ sqrt(pow(shipper[k][2] - customer[j][0], 2) + pow(shipper[k][3] - customer[j][1], 2))\
						- sqrt(pow(shipper[k][2] - shipper[k][0], 2) + pow(shipper[k][3] - shipper[k][1], 2)));
			}
		}
	}
	cout << "\n";

	cout << "cost_y\n";
	for (int j = 0; j < J_customer; j++)
	{
		cost_y[j] = 0.6 * sqrt(pow(customer[j][0], 2) + pow(customer[j][1], 2));
		cout << cost_y[j] << "\t";
	}
	cout << "\n";

	cout << "cost_v\n";
	for (int j = 0; j < J_customer; j++)
	{
		cost_v[j] = 1.5 * cost_y[j];
		cout << cost_v[j] << "\t";
	}
	cout << "\n";

	delete2DArray(mobile, I_mobile);
	delete2DArray(customer, J_customer);
	delete2DArray(shipper, K_shipper);
}

double** Allocate2DArray(int len1, int len2) {
	double** newArray = new double* [len1];
	for (int i = 0; i < len1; i++) {
		newArray[i] = new double[len2];
	}
	return newArray;
}

void delete2DArray(double** array, int len1) {
	for (int i = 0; i < len1; i++) {
		delete[] array[i];
	}
}

void sampling(IloEnv env, int size) {

	sample = IloArray<scenario>(env, size);

	random_device rd;
	mt19937 generator(rd());
	uniform_real_distribution<double> rnum(0.0, 1.0);

	unordered_set<string> subsets; //a temporary set for identifying repeated sample unit

	int sampled = 0;

	while(sampled < size) {
		string sampleUnit = " ";
		scenario tempUnit = scenario(env, K_shipper);
		for (int k = 0; k < K_shipper; k++) {
			double show = rnum(generator);
			if (show < p_cs) {// crowd-shipper shows up
				sampleUnit += "1";
				tempUnit[k] = 1.0;
			}
			else {// crowd-shipper does not show up
				sampleUnit += "0";
				tempUnit[k] = 0.0;
			}
		}
		if (subsets.count(sampleUnit) == 0) {
			subsets.insert(sampleUnit);
			sample[sampled] = scenario(env, K_shipper);
			for (int k = 0; k < K_shipper; k++) sample[sampled][k] = tempUnit[k];
			sampled ++;
		}
		tempUnit.end();
	}
}

void assignSol(IloNumArray z_val, IloNumArray y_val, IloNumArray2D w_val, IloNumArray z_get, IloNumArray y_get, IloNumArray2D w_get) {
	
	for (int i = 0; i < I_mobile; i++) {
		z_get[i] = z_val[i];
		for (int j = 0; j < J_customer; j++) {
			w_get[i][j] = w_val[i][j];
		}
	}
	for (int j = 0; j < J_customer; j++) y_get[j] = y_val[j];
}

void assignSample(IloArray<scenario> sample_get) {
	for (int i = 0; i < sample.getSize(); i++) {
		for (int k = 0; k < K_shipper; k++) {
			sample_get[i][k] = sample[i][k];
		}
	}
}

void MeanAndStd(vector<double> data, double& mean, double& std) {
	mean = accumulate(data.begin(), data.end(), 0.0) / data.size();

	double squaredDiffSum = 0.0;
	for (double value : data) {
		double diff = value - mean;
		squaredDiffSum += diff * diff;
	}
	double variance = squaredDiffSum / (data.size() - 1);
	std = sqrt(variance/ data.size()); //sample standard deviation
}


