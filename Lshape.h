#include "Parameters.h"

#ifndef _Lshape
#define _Lshape
void LShape
(IloEnv env, IloExpr& objFirst, Var1D z, Var1D y, Var2D w, Var1D theta, Var1D v, Var3D p,
	IloNum& obj_opt, IloNumArray z_opt, IloNumArray y_opt, IloNumArray2D w_opt, IloNumArray theta_opt,
	IloExpr& objSecond, IloRangeArray constr9b, IloRangeArray constr9c, Constr2D constr9d);
bool checkBranch(IloCplex cplex, Var1D z, Var1D y, Var2D w, branchNode& branchLeaf);
void fathom(vector<branchNode>& branchTree, Var1D z, Var1D y, Var2D w, bool& intOpt);
void changeBounds(branchNode node, Var1D z, Var1D y, Var2D w, int option);
bool optCut(IloEnv env, IloExpr objSecond, IloRangeArray constr9b, IloRangeArray constr9c, Constr2D constr9d, IloNumVar theta, Var1D y, Var2D w, Var1D v, Var3D p, IloNum theta_val, scenario xi, IloRange& cut);
IloNum secondStageObj(IloEnv env, IloExpr objSecond, IloRangeArray constr9b, IloRangeArray constr9c, Constr2D constr9d, Var1D v, Var3D p, scenario xi);
#endif


