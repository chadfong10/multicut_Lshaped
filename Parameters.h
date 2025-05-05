#include "ilcplex\cplex.h"
#include "ilcplex\ilocplex.h"
#include <vector>
#include <string>
#include <iostream>
#include <tuple>
#include <random>
using namespace std;

#ifndef _Par_H
#define _Par_H
typedef int set_cardinality;
typedef IloNumArray scenario;

typedef struct BRANCH {
	char varTitle;
	int index;
	double value;
	bool secondTraverse = false;
}branchNode;

typedef IloArray<IloNumArray> IloNumArray2D;
typedef IloArray<IloArray<IloNumArray>> IloNumArray3D;

typedef IloNumVarArray Var1D;
typedef IloArray<IloNumVarArray> Var2D;
typedef IloArray<IloArray<IloNumVarArray>> Var3D;

typedef IloIntVarArray IntVar1D;
typedef IloArray<IloIntVarArray> IntVar2D;
typedef IloArray<IloArray<IloIntVarArray>> IntVar3D;
typedef IloArray<IloArray<IloArray<IloIntVarArray>>> IntVar4D;

typedef IloArray<IloRangeArray> Constr2D;
typedef IloArray<Constr2D> Constr3D;

//Sets and parameters
extern bool FindOptimum;
extern set_cardinality I_mobile, J_customer, K_shipper;
extern IloNumArray cost_z;
extern IloNumArray3D cost_p;
extern IloNumArray cost_y;
extern IloNumArray cost_v;
extern IloBoolArray t_mobile, t_shipper;
extern int C_capacity;
extern double p_cs;
extern int sampleSize;
extern int numSamples;
extern int pickSize;
extern int evalSize;
extern IloArray<scenario> sample;
void inputParameters(IloEnv env, string filename);
void sampling(IloEnv env, int size);
void assignSol(IloNumArray z_val, IloNumArray y_val, IloNumArray2D w_val, IloNumArray z_get, IloNumArray y_get, IloNumArray2D w_get);
void assignSample(IloArray<scenario> sample_get);
void MeanAndStd(vector<double> data, double& mean, double& std);
#endif