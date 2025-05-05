#include "Parameters.h"

#ifndef _Msetting
#define _Msetting
void setFirstStageConstr(IloModel modelFirst, Var2D w, Var1D y, Var1D z);
void setSecondStage(IloEnv env, IloExpr & objSecond, IloNumVarArray v, Var3D p, IloRangeArray constr9b, IloRangeArray constr9c, Constr2D constr9d);
#endif
