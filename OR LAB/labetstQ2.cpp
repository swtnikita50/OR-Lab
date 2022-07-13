#include<iostream>
#include<vector>
#include<algorithm>//To use the function next_permutation directly
#include<cmath>//for isnan()

using namespace std;
#define vd vector<double>
#define vvd vector<vector<double>>
#define vi vector<int>
#define num_iterations_max 250
#define num_iterations_min 10
#define pb push_back
#define epsilon_threshold 0.00001
#define double_MAXX 1000000
#define zero_threshold 0.0000001
#define counter_overflow_threshold 15
#define M_value 1000000

//FOR DEBUGGING
////////////////////////////
void pritntvec(vi arr) {
	for (int i = 0; i < arr.size(); i++)
		cout << arr[i] << " ";
	cout << endl;
}

void pritntvecd(vd arr) {
	for (int i = 0; i < arr.size(); i++)
		cout << arr[i] << " ";
	cout << endl;
}

void printmat(vvd mat) {
	for (int i = 0; i < mat.size(); i++) {
		for (int j = 0; j < mat[i].size(); j++)
			cout << mat[i][j] << " ";
		cout << endl;
	}
}
///////////////////////////

double absol(double a) {
	if (a < 0) a *= -1;
	return a;
}

//Find dot product of 2 n dimensional vectors
double find_dot_prod(vd arr1, vd arr2) {
	int n = arr1.size();
	double ret = 0;
	for (int i = 0; i < n; i++)
		ret += (arr1[i] * arr2[i]);
	return ret;
}

//Function to be used to decide when to terminate in solve_gauss_siedel
double find_magnitude(vd arr) {
	double ret = 0;
	for (int i = 0; i < arr.size(); i++) {
		ret += (arr[i] * arr[i]);
	}
	return ret;
}


//Construct a matrix(from the input set of equations) to be used for the Gauss-Siedel method
vvd form_mat(vi index_set, vvd mat) {
	int m = mat.size();
	int n = (mat[0].size()) - 1;
	vvd ret;
	for (int i = 0; i < m; i++) {
		vd temp;
		for (int j = 0; j < n; j++) {
			if (index_set[j] == 0) continue;
			temp.pb(mat[i][j]);
		}
		temp.pb(mat[i][n]);
		ret.pb(temp);
	}

	return ret;
}

//Calculate the value of the objective function for the current values of the variables
double find_obj_fn_value(vd obj_fn, vd ans, vi index_set) {
	int n = index_set.size();
	double ret = 0;
	int count = 0;
	for (int i = 0; i < n; i++) {
		if (index_set[i] == 0)
			continue;

		if (ans[count] < 0 || isnan(ans[count])) {
			//cout<<"This is an infeasible solution"<<endl<<endl;
			return -1;
		}
		ret += (obj_fn[i] * ans[count++]);
	}

	if (ret < 0) return -1;
	count = 0;
	for (int i = 0; i < n; i++) {
		if (index_set[i] == 0) cout << "0(NB_Var) ";//NB_Var = Non-basic variable
		else cout << ans[count++] << " ";
	}

	cout << ",  Objective_fn value: " << ret << endl << endl;
	return ret;
}

//Modify the objective function to take the slack variables into account
void create_maximisation_objective_function(vd& objective_fn, bool is_max, int n, vi slack_list, vi surplus_list, vi artificial_list) {
	for (int i = objective_fn.size(); i < n; i++) objective_fn.pb(0);
	cout << "objective_fn_size: " << objective_fn.size() << endl;
	cout << "n: " << n << endl;
	if (!is_max) {
		for (int i = 0; i < objective_fn.size(); i++) objective_fn[i] = -1 * objective_fn[i];
	}

	//Use artificial variables to modify the objective function to include the bigM terms
	for (int i = 0; i < artificial_list.size(); i++)
		objective_fn[artificial_list[i] - 1] = -1 * M_value;
}

//Create the full coefficient matrix which includes the slack variables as well
void create_mat(vi eqn_types, vvd& mat, vi& slack_list, vi& surplus_list, vi& artificial_list) {
	int m = mat.size();
	int n = mat[0].size();

	vd last_column(m, 0);
	for (int i = 0; i < m; i++) {
		last_column[i] = mat[i][n - 1];
		mat[i].pop_back();
	}
	n--;//n is the number of variables till now

	//Add the slack, surplus, artificial variables to the matrix
	for (int i = 0; i < m; i++) {
		if (eqn_types[i] == 1) { //<=
			slack_list.pb(n + 1); n++;
			for (int j = 0; j < m; j++) {
				if (j == i) mat[j].pb(1);
				else mat[j].pb(0);
			}
		}

		else if (eqn_types[i] == 2) { //>=
			surplus_list.pb(n + 1);
			artificial_list.pb(n + 2);
			n = n + 2;//add surplus and artificial variable
			for (int j = 0; j < m; j++) {
				if (j == i) { mat[j].pb(-1); mat[j].pb(1); }//subtract surplus and add artificial variable
				else { mat[j].pb(0); mat[j].pb(0); }
			}
		}

		else if (eqn_types[i] == 3) { //=
			artificial_list.pb(n + 1); n++;//add artificial variable
			for (int j = 0; j < m; j++) {
				if (j == i) mat[j].pb(1);
				else mat[j].pb(0);
			}
		}
	}

	//Add the RHS of the contraints back again into the matrix
	for (int i = 0; i < m; i++) {
		mat[i].pb(last_column[i]);
	}
}

bool can_simplex_be_used(vi eqn_types) {
	for (int i = 0; i < eqn_types.size(); i++) {
		if (eqn_types[i] != 1) { cout << endl << "Normal Simplex method cant be used to solve this system, using BIG_M method" << endl; return 0; }
	}
	cout << "Using simplex method to solve optimisation question" << endl;
	return 1;
}

vvd construct_simplex_tableau(vvd mat) {
	int num_rows = mat.size();
	int num_cols = mat[0].size();
	vector<vector<double>> tableau(num_rows, vector<double>(num_cols + 2, 0));//2 extra columns for C0, Basis
	for (int i = 0; i < num_rows; i++) {
		tableau[i][0] = 0;
		tableau[i][1] = num_cols - num_rows + i + 1;//column index of the basis variable inside the tableau
	}
	for (int i = 0; i < num_rows; i++) {
		for (int j = 0; j < num_cols; j++) tableau[i][j + 2] = mat[i][j];
	}
	return tableau;
}

//vvd tableau = construct_bigM_tableau(mat, objective_fn, slack_variable_list, surplus_variable_list, artificial_variable_list)
vvd construct_bigM_tableau(vvd mat, vd objective_fn, vi slack_list, vi surplus_list, vi artificial_list) {
	int num_rows = mat.size();
	int num_cols = mat[0].size();
	vector<vector<double>> tableau(num_rows, vector<double>(num_cols + 2, 0));//2 extra columns for C0, Basis

	int curr_row = 0;
	int slack_index = 0;
	int artificial_index = 0;
	while ((slack_index < slack_list.size()) && (artificial_index < artificial_list.size())) {
		if (slack_list[slack_index] < artificial_list[artificial_index]) {
			tableau[curr_row][1] = slack_list[slack_index++] + 1;//column index in tableau
			tableau[curr_row++][0] = 0;
		}
		else {
			tableau[curr_row][1] = artificial_list[artificial_index++] + 1;
			tableau[curr_row++][0] = -1 * M_value;
		}
	}

	while (slack_index < slack_list.size()) {
		tableau[curr_row][1] = slack_list[slack_index++] + 1;//column index in tableau
		tableau[curr_row++][0] = 0;
	}

	while (artificial_index < artificial_list.size()) {
		tableau[curr_row][1] = artificial_list[artificial_index++] + 1;
		tableau[curr_row++][0] = -1 * M_value;
	}


	//construct rest of the tableau
	for (int i = 0; i < num_rows; i++) {
		for (int j = 0; j < num_cols; j++) tableau[i][j + 2] = mat[i][j];
	}
	return tableau;
}


vd calculate_deviations(vvd tableau, vd objective_fn) {
	int num_rows = tableau.size();
	int num_cols = tableau[0].size();
	vd ret(num_cols - 3, 0);

	vd col(num_rows, 0);
	for (int i = 0; i < num_rows; i++) col[i] = tableau[i][0];

	for (int col_index = 2; col_index < num_cols - 1; col_index++) {
		vd temp(num_rows, 0);
		for (int i = 0; i < num_rows; i++) temp[i] = tableau[i][col_index];

		ret[col_index - 2] = objective_fn[col_index - 2] - find_dot_prod(temp, col);
		temp.clear();
	}
	return ret;
}

bool is_termination_reached(vd deviations) {
	for (int i = 0; i < deviations.size(); i++) {
		if (deviations[i] > zero_threshold) return 0;
	}
	return 1;
}

int find_max_elem_index(vd arr) {
	double max = 0;
	int max_index = -1;
	for (int i = 0; i < arr.size(); i++) {
		if (arr[i] > max) {
			max = arr[i];
			max_index = i;
		}
	}
	return max_index;
}

// R1 -> R1 - dR2
void subtract_rows(vvd& tableau, double d, int r1, int r2) {
	int num_rows = tableau.size();
	int num_cols = tableau[0].size();
	for (int i = 2; i < num_cols; i++) tableau[r1][i] = tableau[r1][i] - d * tableau[r2][i];
}

// R1 -> R1/d
void divide_row(vvd& tableau, double d, int r) {
	int num_rows = tableau.size();
	int num_cols = tableau[0].size();
	for (int i = 2; i < num_cols; i++) tableau[r][i] = tableau[r][i] / d;
}

bool do_multiple_solutions_exist(vd deviations, int num_rows) {
	int num_zeros = 0;
	for (int i = 0; i < deviations.size(); i++) {
		if ((deviations[i] > (-1 * zero_threshold)) && (deviations[i] < zero_threshold)) num_zeros++;
	}
	if (num_zeros > num_rows) return 1;
	return 0;
}

bool is_in(int d, vi arr) {
	for (int i = 0; i < arr.size(); i++) { if (arr[i] == d) return 1; }
	return 0;
}

void solve_problem(vvd tableau, vd objective_fn, bool is_max, vi artificial_list) {
	int num_rows = tableau.size();
	int num_cols = tableau[0].size();

	cout << endl; int count = 0;
	while (1) {
		if (count > counter_overflow_threshold) { cout << "Excess iterations, returning!!...." << endl; return; }
		vd deviations = calculate_deviations(tableau, objective_fn);
		if (count == 1) {
			cout << endl << endl << "Iteration: " << count + 1 << endl;
			//cout << "The Tableau: " << endl; printmat(tableau);


			//cout << "Zj-Cj: "; pritntvecd(deviations);
		}
		bool is_end = is_termination_reached(deviations);
		if (is_end) { cout << "Reached the termination state" << endl; break; }

		//Identifying the entering variable
		int entering_var_col_index = find_max_elem_index(deviations) + 2;
		if (entering_var_col_index == -1) { cout << "Something is Wrong!" << endl; break; }

		//Identifying the leaving variable
		int leaving_var_row_index = -1;
		double min_ratio = double_MAXX;
		for (int i = 0; i < num_rows; i++) {
			if (tableau[i][entering_var_col_index] <= 0) continue;
			double temp = tableau[i][num_cols - 1] / tableau[i][entering_var_col_index];
			if (temp < min_ratio) {
				min_ratio = temp;
				leaving_var_row_index = i;
			}
		}
		if (leaving_var_row_index == -1) { cout << "Unbounded solution!!!....EXITING..." << endl; return; }
		if (count == 1) {
			//cout << "leaving_var_row_index: " << leaving_var_row_index << ", entering_var_col_index: " << entering_var_col_index << endl;
			cout << "The pivot row elements and their values are: " << endl;
			for (int l = 2; l < tableau[0].size(); l++) {
				cout << "x[" << l-1 << "]= " << tableau[leaving_var_row_index][l] << endl;
				}
			cout << endl << endl;
		}

		//Perform row operations on the tableau to get the tableau for the next iteration ready
		for (int row_index = 0; row_index < num_rows; row_index++) {
			if (row_index == leaving_var_row_index) divide_row(tableau, tableau[row_index][entering_var_col_index], row_index);
			else {
				double ratio_temp = tableau[row_index][entering_var_col_index] / tableau[leaving_var_row_index][entering_var_col_index];
				subtract_rows(tableau, ratio_temp, row_index, leaving_var_row_index);
			}
		}

		//Modify 1st and 2nd columns of the tableau
		tableau[leaving_var_row_index][0] = objective_fn[entering_var_col_index - 2];
		tableau[leaving_var_row_index][1] = entering_var_col_index;

		deviations.clear();
		count++;
	}

	int num_vars = num_cols - 3;
	vd variable_vals(num_vars, 0);
	for (int i = 0; i < num_rows; i++)	variable_vals[tableau[i][1] - 2] = tableau[i][num_cols - 1];

	//Identify infeasible soution
	for (int i = 0; i < num_rows; i++) {
		if (is_in(tableau[i][1] - 1, artificial_list)) {
			cout << endl << "The iterations have been completed and there are artificial variables in the base with values strictly greater than 0, so the problem has no solution (INFEASIBLE!!!!!)." << endl;
			return;
		}
	}

	//cout << endl << "The Optimisation function is optimised at the point : ";
	for (int i = 0; i < num_rows; i++) {
		//cout << "x" << tableau[i][1] - 1 << " = " << tableau[i][num_cols - 1] << ", ";
	}
	//cout << "Rest all xis=0" << endl;

	double optimal_value = find_dot_prod(objective_fn, variable_vals);
	if (!is_max) optimal_value *= -1;
	//cout << "The value of the objective function at this point is : " << optimal_value << endl;


	//CHECKING FOR MULTIPLE SOLUTIONS
	vd deviations = calculate_deviations(tableau, objective_fn);
	if (do_multiple_solutions_exist(deviations, tableau.size())) {
		cout << "We are at an optimal point and there are non-basic variables with reduced cost equal to 0, so there are multiple values for the decision variables that allow obtaining the optimal value" << endl;
		cout << "Along with the mentioned point, infinite other points exist s.t. the optimisation function is optimised" << endl;
		return;
	}
}


int main() {
	int n, m;

	cout << "No. of Constraints(m): "; cin >> m;//num_eqns
	cout << "No. of Variables(n): "; cin >> n;//num_vars
	vector<vector<double>> mat(m, vector<double>(n + 1, 0));

	for (int i = 0; i < m; i++) {
		cout << "Enter Constraint " << i << " coefficients: " << endl;
		for (int j = 0; j < n + 1; j++) {
			cin >> mat[i][j];
		}
	}
	cout << endl << "Enter whether the equations are <= type(1), >= type(2), or = type(3): " << endl;
	vi eqn_types(m, 0);
	for (int i = 0; i < m; i++) {
		cout << "Equation " << i << ": "; cin >> eqn_types[i];
	}
	bool is_problem_simplex = can_simplex_be_used(eqn_types);
	vi slack_variable_list; vi surplus_variable_list; vi artificial_variable_list;
	create_mat(eqn_types, mat, slack_variable_list, surplus_variable_list, artificial_variable_list);

	printmat(mat);
	cout << "slack variables: "; pritntvec(slack_variable_list);
	cout << "surplus variables: "; pritntvec(surplus_variable_list);
	cout << "artifiical variables: "; pritntvec(artificial_variable_list);

	cout << endl << "Enter the objective function: " << endl;
	vd objective_fn(n, 0);
	for (int i = 0; i < n; i++) {
		cin >> objective_fn[i];
	}
	bool is_maximisation;
	cout << endl << "Does the objective_fn have to be maximised(1) or minimised(0)? : ";
	cin >> is_maximisation;
	n = mat[0].size();
	create_maximisation_objective_function(objective_fn, is_maximisation, n - 1, slack_variable_list, surplus_variable_list, artificial_variable_list);

	//cout << "The coeficient matrix: " << endl;
	//printmat(mat);
	//cout << "The objective function: " << endl;
	//pritntvecd(objective_fn);
	cout << endl;

	cout << "Solving Using BigM Method..." << endl << "Value of M is: " << M_value << endl;

	if (is_problem_simplex) {
		vvd tableau = construct_simplex_tableau(mat);
		//cout << "The Simplex tableau: " << endl;
		//printmat(tableau);
		solve_problem(tableau, objective_fn, is_maximisation, artificial_variable_list);
	}


	else {
		vvd tableau = construct_bigM_tableau(mat, objective_fn, slack_variable_list, surplus_variable_list, artificial_variable_list);
		//cout << "The bigM tableau: " << endl;
		//printmat(tableau);
		solve_problem(tableau, objective_fn, is_maximisation, artificial_variable_list);
	}
}

/*
Input for Question:
2
4
5 4 2 1 100
2 3 8 1 75
2
2
12 8 14 10
0
*/
