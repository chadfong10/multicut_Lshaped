#include "Lshape.h"
#include "modelSetting.h"

void LShape
(IloEnv env, IloExpr& objFirst, Var1D z, Var1D y, Var2D w, Var1D theta, Var1D v, Var3D p,
	IloNum& obj_opt, IloNumArray z_opt, IloNumArray y_opt, IloNumArray2D w_opt, IloNumArray theta_opt,
	IloExpr& objSecond, IloRangeArray constr9b, IloRangeArray constr9c, Constr2D constr9d) {

	IloModel modelFirst(env);
	IloCplex cplexFirst(modelFirst);
	cplexFirst.setOut(env.getNullStream());
	cplexFirst.setParam(IloCplex::TiLim, 1.0);
	cplexFirst.setParam(IloCplex::EpAGap, 1.0e-6);

	modelFirst.add(IloMinimize(env, objFirst));

	//set up the basic constraints of the first-stage model
	setFirstStageConstr(modelFirst, w, y, z);

	/*L-shape*/
	IloNum incumbent = IloInfinity;
	bool intOptFirst = false;
	vector<branchNode> branchTree;
	bool optimum = false;
	IloRange cut;

	int iter = 0;

	while (!optimum) {

		//solve the first stage
		intOptFirst = false;
		incumbent = IloInfinity;
		while (!intOptFirst) {
			cplexFirst.solve();
			if (cplexFirst.getStatus() == IloAlgorithm::Optimal) {
				branchNode branchLeaf;
				if (cplexFirst.getObjValue() < incumbent) {
					if (checkBranch(cplexFirst, z, y, w, branchLeaf)) {
						//branching
						branchTree.push_back(branchLeaf);
						changeBounds(branchLeaf, z, y, w, 1);
					}
					else {
						//New incumbent, and then we move backward on the tree(fathom)
						incumbent = cplexFirst.getObjValue();
						cplexFirst.getValues(z, z_opt);
						cplexFirst.getValues(y, y_opt);
						for (int i = 0; i < I_mobile; i++) {
							cplexFirst.getValues(w[i], w_opt[i]);
						}
						cplexFirst.getValues(theta, theta_opt);
						fathom(branchTree, z, y, w, intOptFirst);
					}
				}
				else {//Larger than incumbent so fathom
					fathom(branchTree, z, y, w, intOptFirst);
				}
			}
			else if (cplexFirst.getStatus() == IloAlgorithm::Infeasible) {
				fathom(branchTree, z, y, w, intOptFirst);
			}
			else {
				cout << "Unbounded status or other unknown CPLEX status appears.\n";
				exit(1);
			}
		}
		// release the memory space
		branchTree.shrink_to_fit();

		/*Optimality cut*/
		//change the right-hand side of constraints in second stage
		for (int j = 0; j < J_customer; j++) {
			constr9c[j].setBounds(1.0 - y_opt[j], 1.0 - y_opt[j]);
		}
		for (int i = 0; i < I_mobile; i++) {
			constr9d[i].setUbs(w_opt[i]);
		}

		//Optimality cut
		optimum = true;
		for (int s = 0; s < sampleSize; s++) {
			if (optCut(env, objSecond, constr9b, constr9c, constr9d, theta[s], y, w, v, p, theta_opt[s], sample[s], cut)) {
				modelFirst.add(cut);
				optimum = false;
			}
		}
		if (!(iter % 10)) {
			cout << "Objective value in iteration " << iter << " : " << incumbent << endl;
		}
		iter++;
	}
	obj_opt = incumbent;
	cplexFirst.end();
	modelFirst.end();
}

//Check integer constraints satisfication and apply branch and bound
bool checkBranch(IloCplex cplex, Var1D z, Var1D y, Var2D w, branchNode& branchLeaf) {
	IloEnv env = cplex.getEnv();

	bool whetherBranch = false;

	//consider z
	IloNumArray z_result = IloNumArray(env, I_mobile);
	cplex.getValues(z_result, z);
	for (int i = 0; i < I_mobile; i++) {
		if (z_result[i] < 1.0 && z_result[i] > 0.0) {
			whetherBranch = true;
			branchLeaf.varTitle = 'z';
			branchLeaf.index = i;
			branchLeaf.value = 0.0;
			break;
		}
	}

	if (!whetherBranch) {
		//consider y
		IloNumArray y_result = IloNumArray(env, J_customer);
		cplex.getValues(y_result, y);
		for (int j = 0; j < J_customer; j++) {
			if (y_result[j] < 1.0 && y_result[j] > 0.0) {
				whetherBranch = true;
				branchLeaf.varTitle = 'y';
				branchLeaf.index = j;
				branchLeaf.value = 1.0;
				break;
			}
		}
	}

	if (!whetherBranch) {
		//consider w
		for (int i = 0; i < I_mobile; i++) {
			IloNumArray w_result = IloNumArray(env, J_customer);
			cplex.getValues(w_result, w[i]);
			for (int j = 0; j < J_customer; j++) {
				if (w_result[j] < 1.0 && w_result[j] > 0.0) {
					whetherBranch = true;
					branchLeaf.varTitle = 'w';
					branchLeaf.index = i * J_customer + j;
					branchLeaf.value = 0.0;
					break;
				}
			}
		}
	}

	return whetherBranch;
}

void fathom(vector<branchNode>& branchTree, Var1D z, Var1D y, Var2D w, bool& intOpt) {
	if (branchTree.empty()) { // If the first solution is already an integer-feasible solution.
		intOpt = true;
	}
	else {
		branchNode preBranchNode = branchTree.back();
		branchTree.pop_back();
		while (!branchTree.empty() && preBranchNode.secondTraverse) {
			changeBounds(preBranchNode, z, y, w, 3); // after second traverse
			preBranchNode = branchTree.back();
			branchTree.pop_back();
		}
		if (!preBranchNode.secondTraverse) { // second traverse
			branchTree.size();
			preBranchNode.secondTraverse = true;
			branchTree.push_back(preBranchNode);
			changeBounds(preBranchNode, z, y, w, 2);
		}
		else {
			//We clear the stack, so we have found the optimum
			changeBounds(preBranchNode, z, y, w, 3);
			intOpt = true;
		}
	}
}

void changeBounds(branchNode node, Var1D z, Var1D y, Var2D w, int option) {
	switch(option){
		case (1): { //new branch
			if (node.varTitle == 'z') {
				z[node.index].setBounds(node.value, node.value);
			}
			else if (node.varTitle == 'y') {
				y[node.index].setBounds(node.value, node.value);
			}
			else { //w
				w[node.index / J_customer][node.index % J_customer].setBounds(node.value, node.value);
			}
			break;
		}
		case (2): { //previously branched (second traverse)
			if (node.varTitle == 'z') {
				z[node.index].setBounds(1.0 - node.value, 1.0 - node.value);
			}
			else if (node.varTitle == 'y') {
				y[node.index].setBounds(1.0 - node.value, 1.0 - node.value);
			}
			else { //w
				w[node.index / J_customer][node.index % J_customer].setBounds(1.0 - node.value, 1.0 - node.value);
			}
			break;
		}
		case (3): { //reset the bounds to (0.0, 1.0) after second traverse
			if (node.varTitle == 'z') {
				z[node.index].setBounds(0.0, 1.0);
			}
			else if (node.varTitle == 'y') {
				y[node.index].setBounds(0.0, 1.0);
			}
			else { //w
				w[node.index / J_customer][node.index % J_customer].setBounds(0.0, 1.0);
			}
			break;
		}
	}
}

bool optCut(IloEnv env, IloExpr objSecond, IloRangeArray constr9b, IloRangeArray constr9c, Constr2D constr9d, IloNumVar theta, Var1D y, Var2D w, Var1D v, Var3D p, IloNum theta_val, scenario xi, IloRange& cut) {
	
	//change the right-hand sides of constraints 9b in second stage
	constr9b.setUbs(xi);

	// Second stage model of a scenario
	IloModel modelSecond(env);
	IloCplex cplexSecond(modelSecond);
	cplexSecond.setOut(env.getNullStream());
	cplexSecond.setParam(IloCplex::TiLim, 1.0);
	cplexSecond.setParam(IloCplex::EpAGap, 1.0e-6);


	//objective
	modelSecond.add(IloMinimize(env, objSecond));

	//constraints 9b
	modelSecond.add(constr9b);

	//constraint 9c
	modelSecond.add(constr9c);

	//constraint 9d
	for (int i = 0; i < I_mobile; i++) {
		modelSecond.add(constr9d[i]);
	}

	cplexSecond.solve();

	//Check the totally unimodular property
	/*IloNumArray v_val = IloNumArray(env, J_customer);
	cplexSecond.getValues(v, v_val);
	for (int j = 0; j < J_customer; j++) {
		if (v_val[j] > 0.0 && v_val[j] < 1.0) {
			cout << "Not totally unimodular!!!!!\n";
			cout << "v_" << j << ":" << v_val[j];
			exit(1);
		}
	}
	IloNumArray p_val = IloNumArray(env, K_shipper);
	for (int i = 0; i < I_mobile; i++) {
		for (int j = 0; j < J_customer; j++) {
			cplexSecond.getValues(p[i][j], p_val);
			for (int k = 0; k < K_shipper; k++) {
				if (p_val[k] > 0.0 && p_val[k] < 1.0) {
					cout << "Not totally unimodular!!!!!\n";
					cout << "p_" << i << "_" << j << "_" << k << ":" << p_val[k];
					exit(1);
				}
			}
		}
	}*/


	//compare objective and theta
	if (theta_val >= cplexSecond.getObjValue() - 1.0e-4) {
		cplexSecond.end();
		modelSecond.end();
		return false;
	}

	/*Dual*/
	IloNumArray dual9b(env, K_shipper), dual9c(env, J_customer);
	cplexSecond.getDuals(dual9b, constr9b);
	cplexSecond.getDuals(dual9c, constr9c);
	IloArray<IloNumArray> dual9d(env, I_mobile);
	for (int i = 0; i < I_mobile; i++) {
		dual9d[i] = IloNumArray(env, J_customer);
		cplexSecond.getDuals(dual9d[i], constr9d[i]);
	}

	/*cut*/
	IloExpr cutExpr =  - theta + IloScalProd(dual9b, xi) + IloSum(dual9c) - IloScalProd(dual9c, y);
	for (int i = 0; i < I_mobile; i++) {
		cutExpr += IloScalProd(dual9d[i], w[i]);
	}
	cut = IloRange(cutExpr <= 0);

	cplexSecond.end();
	modelSecond.end();
	return true;
}

IloNum secondStageObj(IloEnv env, IloExpr objSecond, IloRangeArray constr9b, IloRangeArray constr9c, Constr2D constr9d, Var1D v, Var3D p, scenario xi) {

	//change the right-hand side of constraints in second stage
	constr9b.setUbs(xi);

	// Second stage model of a scenario
	IloModel modelSecond(env);
	IloCplex cplexSecond(modelSecond);
	cplexSecond.setOut(env.getNullStream());
	cplexSecond.setParam(IloCplex::EpAGap, 1.0e-4);

	//objective
	modelSecond.add(IloMinimize(env, objSecond));

	//constraints 9b
	modelSecond.add(constr9b);

	//constraint 9c
	modelSecond.add(constr9c);

	//constraint 9d
	for (int i = 0; i < I_mobile; i++) {
		modelSecond.add(constr9d[i]);
	}

	cplexSecond.solve();

	IloNum obj_opt = cplexSecond.getObjValue();

	cplexSecond.end();
	modelSecond.end();
	return obj_opt;
}