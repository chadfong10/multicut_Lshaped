#include "modelSetting.h"

void
setFirstStageConstr(IloModel modelFirst, Var2D w, Var1D y, Var1D z) {
	IloEnv env = modelFirst.getEnv();

	//constraints
	//7b
	/*for (int i = 0; i < I_mobile; i++) {
		for (int j = 0; j < J_customer; j++) {
			modelFirst.add(IloRange(w[i][j] - z[i] <= 0));
		}
	}*/

	//7c
	for (int j = 0; j < J_customer; j++) {
		IloExpr sum(env);
		for (int i = 0; i < I_mobile; i++) {
			sum += w[i][j];
		}
		modelFirst.add(IloRange(sum + y[j] == 1));
		sum.end();
	}

	//7d
	/*for (int i = 0; i < I_mobile; i++) {
		modelFirst.add(IloRange(IloSum(w[i]) <= C_capacity));
	}*/

	//combine 7b and 7d
	for (int i = 0; i < I_mobile; i++) {
		modelFirst.add(IloRange(IloSum(w[i]) - C_capacity * z[i] <= 0));
	}

	//valid inequality
	modelFirst.add(IloRange(C_capacity * IloSum(z) - (J_customer - IloSum(y)) >= 0));
}

void setSecondStage(IloEnv env, IloExpr& objSecond, IloNumVarArray v, Var3D p, IloRangeArray constr9b, IloRangeArray constr9c, Constr2D constr9d) {
	//objective
	objSecond = IloScalProd(cost_v, v);
	for (int i = 0; i < I_mobile; i++) {
		for (int j = 0; j < J_customer; j++) {
			objSecond += IloScalProd(cost_p[i][j], p[i][j]);
		}
	}

	//constraints 9b
	for (int k = 0; k < K_shipper; k++) {
		IloExpr sum(env);
		for (int i = 0; i < I_mobile; i++) {
			for (int j = 0; j < J_customer; j++) {
				sum += p[i][j][k];
			}
		}
		constr9b[k] = IloRange(sum <= 0.0);
		sum.end();
	}

	//constraint 9c
	for (int j = 0; j < J_customer; j++) {
		IloExpr sum(env);
		for (int i = 0; i < I_mobile; i++) sum += IloSum(p[i][j]);
		constr9c[j] = IloRange(sum + v[j] == 1.0);
		sum.end();
	}

	//constraint 9d
	for (int i = 0; i < I_mobile; i++) {
		constr9d[i] = IloRangeArray(env, J_customer);
		for (int j = 0; j < J_customer; j++) {
			constr9d[i][j] = IloRange(IloSum(p[i][j]) <= 0.0);
		}
	}
}