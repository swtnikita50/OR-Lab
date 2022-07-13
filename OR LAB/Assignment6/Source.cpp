//Assignment6 submission
//Use the DUAL SIMPLEX method to solve a linear optimisation problem with constraints of the form <=
//Note: For your convenience you can copy paste the question inputs from the end of the program to test the code

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

bool is_infeasible = 0;

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

//If RHS of any constraint is negative multpily the constraint by -1
void convert_to_standard_form(vvd& mat, vi& eqn_types) {
	int n = mat[0].size() - 1;
	for (int i = 0; i < mat.size(); i++) {
		if (mat[i][n] < 0) {
			for (int j = 0; j < n + 1; j++) mat[i][j] *= -1;
			if (eqn_types[i] == 1) eqn_types[i] = 2;
			else if (eqn_types[i] == 2) eqn_types[i] = 1;
		}
	}
}

bool can_simplex_be_used(vi eqn_types) {
	for (int i = 0; i < eqn_types.size(); i++) {
		if (eqn_types[i] != 1) { cout << endl << "Normal Simplex method cant be used to solve this system, using 2 Phase Simplex method " << endl; return 0; }
	}
	cout << "Using simplex method to solve optimisation question" << endl;
	return 1;
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

//Calculate the value of the cost function for the current values of the variables
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

	cout << ",  cost_fn value: " << ret << endl << endl;
	return ret;
}

//Modify the cost function to take the slack variables into account
void create_maximisation_cost_function(vd& cost_fn, bool is_max, int n) {
	for (int i = cost_fn.size(); i < n; i++) cost_fn.pb(0);
	if (is_max) return;
	cout << "cost_fn_size: " << cost_fn.size() << endl;
	cout << "n: " << n << endl;
	for (int i = 0; i < cost_fn.size(); i++) cost_fn[i] = -1 * cost_fn[i];
}

//Create the full coefficient matrix which includes the slack variables as well
void create_mat(vi eqn_types, vvd& mat) {
	int m = mat.size();
	int n = mat[0].size();

	vd last_column(m, 0);
	for (int i = 0; i < m; i++) {
		last_column[i] = mat[i][n - 1];
		mat[i].pop_back();
	}

	//Add the slack variables to the matrix
	for (int i = 0; i < m; i++) {
		if (eqn_types[i] == 1) {
			n++;
			for (int j = 0; j < m; j++) {
				if (j == i) mat[j].pb(1);
				else mat[j].pb(0);
			}
		}

		else if (eqn_types[i] == 2) {
			n++;
			for (int j = 0; j < m; j++) {
				if (j == i) mat[j].pb(-1);
				else mat[j].pb(0);
			}
		}
	}

	for (int i = 0; i < m; i++) {
		mat[i].pb(last_column[i]);
	}

	for (int i = 0; i < m; i++) {
		if (eqn_types[i] == 2) {
			for (int j = 0; j < mat[i].size(); j++)
				mat[i][j] *= -1;
		}
		eqn_types[i] = 1;
	}
}

vvd construct_tableau(vvd mat) {
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

vd calculate_deviations(vvd tableau, vd cost_fn) {
	int num_rows = tableau.size();
	int num_cols = tableau[0].size();
	vd ret(num_cols - 3, 0);

	vd col(num_rows, 0);
	for (int i = 0; i < num_rows; i++) col[i] = tableau[i][0];

	for (int col_index = 2; col_index < num_cols - 1; col_index++) {
		vd temp(num_rows, 0);
		for (int i = 0; i < num_rows; i++) temp[i] = tableau[i][col_index];

		ret[col_index - 2] = cost_fn[col_index - 2] - find_dot_prod(temp, col);
		temp.clear();
	}
	return ret;
}

//If all RHS entries are positive then we have reached the termination
bool is_termination_reached(vvd tableau) {
	int m = tableau.size();
	int n = tableau[0].size();
	for (int i = 0; i < m; i++) {
		if (tableau[i][n - 1] < (-1 * zero_threshold)) return 0;
	}
	return 1;
}

//If all the deviations are negative
bool is_termination_simplex_reached(vd deviations) {
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


void solve_with_simplex(vvd tableau, vd cost_fn, bool is_maximisation) {
	int num_rows = tableau.size();
	int num_cols = tableau[0].size();

	cout << endl << endl << "Solvable by normal Simplex method... Beginning Simplex method!!" << endl; int count = 0;
	while (1) {
		if (count > counter_overflow_threshold) { cout << "Excess iterations, returning!!...." << endl; return; }
		cout << endl << endl << "Iteration: " << count << endl;
		cout << "The Tableau: " << endl; printmat(tableau);

		vd deviations = calculate_deviations(tableau, cost_fn);
		cout << "Deviations: "; pritntvecd(deviations);
		bool is_end = is_termination_simplex_reached(deviations);
		if (is_end) { cout << "Reached the termination state since all deviations are non positive" << endl; break; }

		//Identifying the incoming variable
		int incoming_var_col_index = find_max_elem_index(deviations) + 2;
		if (incoming_var_col_index == -1) { cout << "LAPHDA!!" << endl; break; }

		//Identifying the outgoing variable
		int outgoing_var_row_index = -1;
		double min_ratio = double_MAXX;
		for (int i = 0; i < num_rows; i++) {
			if (tableau[i][incoming_var_col_index] <= 0) continue;
			double temp = tableau[i][num_cols - 1] / tableau[i][incoming_var_col_index];
			if (temp < min_ratio) {
				min_ratio = temp;
				outgoing_var_row_index = i;
			}
		}
		if (outgoing_var_row_index == -1) { cout << "Unbounded solution!!!....EXITING..." << endl; return; }
		cout << "outgoing_var_row_index: " << outgoing_var_row_index << ", incoming_var_col_index: " << incoming_var_col_index << endl;


		//Perform row operations on the tableau to get the tableau for the next iteration ready
		for (int row_index = 0; row_index < num_rows; row_index++) {
			if (row_index == outgoing_var_row_index) divide_row(tableau, tableau[row_index][incoming_var_col_index], row_index);
			else {
				double ratio_temp = tableau[row_index][incoming_var_col_index] / tableau[outgoing_var_row_index][incoming_var_col_index];
				subtract_rows(tableau, ratio_temp, row_index, outgoing_var_row_index);
			}
		}

		//Modify 1st and 2nd columns of the tableau
		tableau[outgoing_var_row_index][0] = cost_fn[incoming_var_col_index - 2];
		tableau[outgoing_var_row_index][1] = incoming_var_col_index;

		deviations.clear();
		count++;
	}

	int num_vars = num_cols - 3;
	vd variable_vals(num_vars, 0);
	for (int i = 0; i < num_rows; i++)	variable_vals[tableau[i][1] - 2] = tableau[i][num_cols - 1];

	cout << endl << "The Optimisation function is optimised at the point : ";
	for (int i = 0; i < num_rows; i++) {
		cout << "x" << tableau[i][1] - 1 << " = " << tableau[i][num_cols - 1] << ", ";
	}
	cout << "Rest all xis=0" << endl;

	double optimal_value = find_dot_prod(cost_fn, variable_vals);
	if (!is_maximisation) optimal_value *= -1;
	cout << "The value of the cost function at this point is : " << optimal_value << endl;

	//CHECKING FOR MULTIPLE SOLUTIONS
	vd deviations = calculate_deviations(tableau, cost_fn);
	if (do_multiple_solutions_exist(deviations, tableau.size())) {
		cout << "ALTERNATIVE OPTIMUM EXISTS......We are at an optimal point and there are non-basic variables with reduced cost equal to 0, so there are multiple values for the decision variables that allow obtaining the optimal value" << endl;
		cout << "Along with the mentioned point, infinite other points exist s.t. the optimisation function is optimised" << endl;
		return;
	}
}

void solve_with_dual_simplex(vvd tableau, vd cost_fn, bool is_maximisation) {
	int num_rows = tableau.size();
	int num_cols = tableau[0].size();

	cout << endl << endl << "Beginning Dual Simplex method!!" << endl; int count = 0;
	while (1) {
		if (count > counter_overflow_threshold) { cout << "Excess iterations, returning!!...." << endl; return; }
		cout << endl << endl << "Iteration: " << count << endl;
		cout << "The Tableau: " << endl; printmat(tableau);

		bool is_end = is_termination_reached(tableau);
		if (is_end) { cout << "Reached the termination state since all RHS values are positive" << endl; break; }

		vd deviations = calculate_deviations(tableau, cost_fn);
		cout << "Deviations: "; pritntvecd(deviations);

		//Identifying the outgoing variable
		int outgoing_var_row_index = -1;
		double min_till_now = zero_threshold;
		for (int i = 0; i < num_rows; i++) {
			if (tableau[i][num_cols - 1] < min_till_now) {
				min_till_now = tableau[i][num_cols - 1];
				outgoing_var_row_index = i;
			}
		}
		if (outgoing_var_row_index == -1) { cout << "UNBOUNDED!! {Since we could not find a valid outgoing variable}" << endl; return; }
		cout << "Outgoing variable: x" << tableau[outgoing_var_row_index][1] - 1 << endl;

		//Identifying the incoming variable
		double min_positive_ratio_till_now = double_MAXX;
		int incoming_var_col_index = -1;
		for (int i = 2; i < num_cols - 1; i++) {
			if (tableau[outgoing_var_row_index][i] < (-1 * zero_threshold)) {
				double temp_ratio = deviations[i - 2] / tableau[outgoing_var_row_index][i];
				if (temp_ratio > 0 && temp_ratio < min_positive_ratio_till_now) {
					min_positive_ratio_till_now = temp_ratio;
					incoming_var_col_index = i;
				}
			}
		}
		if (incoming_var_col_index == -1) {
			cout << "INFEASIBLE!! {Since we could not find valid incoming variable}" << endl;
			is_infeasible = 1;
			return;
		}
		cout << "Incoming variable: x" << incoming_var_col_index - 1 << endl;

		int num_var = num_cols - 3;
		vd variable_val(num_var, 0);
		for (int i = 0; i < num_rows; i++)	variable_val[tableau[i][1] - 2] = tableau[i][num_cols - 1];

		double optimization_value = find_dot_prod(cost_fn, variable_val);
		cout << "The Value of the Optimization Fucntion is: " << optimization_value << endl;


		//Perform row operations on the tableau to get the tableau for the next iteration ready
		for (int row_index = 0; row_index < num_rows; row_index++) {
			if (row_index == outgoing_var_row_index) divide_row(tableau, tableau[row_index][incoming_var_col_index], row_index);
			else {
				double ratio_temp = tableau[row_index][incoming_var_col_index] / tableau[outgoing_var_row_index][incoming_var_col_index];
				subtract_rows(tableau, ratio_temp, row_index, outgoing_var_row_index);
			}
		}

		//Modify 1st and 2nd columns of the tableau
		tableau[outgoing_var_row_index][0] = cost_fn[incoming_var_col_index - 2];
		tableau[outgoing_var_row_index][1] = incoming_var_col_index;

		deviations.clear();
		count++;
	}

	int num_vars = num_cols - 3;
	vd variable_vals(num_vars, 0);
	for (int i = 0; i < num_rows; i++)	variable_vals[tableau[i][1] - 2] = tableau[i][num_cols - 1];

	cout << endl << "The Optimisation function is optimised at the point : ";
	for (int i = 0; i < num_rows; i++) {
		cout << "x" << tableau[i][1] - 1 << " = " << tableau[i][num_cols - 1] << ", ";
	}
	cout << "Rest all xis=0" << endl;

	double optimal_value = find_dot_prod(cost_fn, variable_vals);

	if (!is_maximisation) optimal_value *= -1;
	cout << "The value of the cost function at this point is : " << optimal_value << endl;


	//CHECKING FOR MULTIPLE SOLUTIONS
	vd deviations = calculate_deviations(tableau, cost_fn);
	if (do_multiple_solutions_exist(deviations, tableau.size())) {
		cout << "ALTERNATIVE OPTIMUM EXISTS......We are at an optimal point and there are non-basic variables with reduced cost equal to 0, so there are multiple values for the decision variables that allow obtaining the optimal value" << endl;
		cout << "Along with the mentioned point, infinite other points exist s.t. the optimisation function is optimised" << endl;
		return;
	}
}

vvd construct_dual_tableau(vvd mat, vd& cost_fn, bool& is_maximisation, vi& eqn_types) {
	int num_rows = mat.size();
	int num_cols = mat[0].size();
	vector<vector<double>> ret_mat(num_cols - 1, vector<double>(num_rows + 1, 0));

	for (int i = 0; i < num_rows; i++) {
		for (int j = 0; j < num_cols - 1; j++) {
			ret_mat[j][i] = mat[i][j];
		}
	}
	for (int i = 0; i < cost_fn.size(); i++) {
		ret_mat[i][num_rows] = cost_fn[i];
	}

	cost_fn.clear();
	for (int i = 0; i < num_rows; i++) {
		cost_fn.push_back(mat[i][num_cols - 1]);
	}

	for (int i = 0; i < eqn_types.size(); i++) {
		if (eqn_types[i] == 1) eqn_types[i] = 2;
		else if (eqn_types[i] == 2) eqn_types[i] = 1;
	}

	is_maximisation = (!is_maximisation);
	return ret_mat;
}

vvd new_mat(vvd mat) {
	vvd ret_mat(mat.size(), vd(mat[0].size(), 0));
	for (int i = 0; i < mat.size(); i++) {
		for (int j = 0; j < mat[0].size(); j++) {
			ret_mat[i][j] = mat[i][j];
		}
	}
	return ret_mat;
}

void make_vec(vd& arr, int n) {
	while (arr.size() > n) arr.pop_back();
}

int main() {
	int n, m;
	cout << "Enter no. of constraints(m): "; cin >> m;
	cout << "Enter no. of variables(n): "; cin >> n;
	int num_vars_og = n;
	vvd raw_matrix(m, vd(n + 1, 0));

	for (int i = 0; i < m; i++) {
		cout << "Enter Constraint " << i << " coefficients: " << endl;
		for (int j = 0; j < n + 1; j++) {
			cin >> raw_matrix[i][j];
		}
	}
	cout << endl << "Enter Equation types: <= type(1), >= type(2), or = type(3): " << endl;
	vi eqn_types(m, 0);
	for (int i = 0; i < m; i++) {
		cout << "Equation " << i << ": "; cin >> eqn_types[i];
	}
	vi var_types(n, 0);
	cout << endl << "Enter the feasible region (for >=0 type(0), R(real no.) type(1), <=0 type(2) :" << endl;
	for (int i = 0; i < n; i++) {
		cout << "Variable x" << i << ": "; cin >> var_types[i];
	}



	//Handle the cases when the RHS of an equation is negaitve
	convert_to_standard_form(raw_matrix, eqn_types);
	bool is_problem_simplex = can_simplex_be_used(eqn_types);

	vvd matrix = new_mat(raw_matrix);
	create_mat(eqn_types, matrix);

	cout << endl << "Enter the cost function: " << endl;
	vd cost_fn(n, 0);
	for (int i = 0; i < n; i++) {
		cin >> cost_fn[i];
	}
	bool is_maximisation;
	cout << endl << "Does the cost_fn have to be maximised(1) or minimised(0)? : ";
	cin >> is_maximisation;
	n = matrix[0].size();
	create_maximisation_cost_function(cost_fn, is_maximisation, n - 1);

	cout << "The coeeficient matrix: " << endl;
	printmat(matrix);
	cout << "The cost function: " << endl;
	pritntvecd(cost_fn);
	cout << endl;

	vvd tableau = construct_tableau(matrix);
	cout << "The Simplex tableau: " << endl;
	printmat(tableau);

	if (is_problem_simplex) solve_with_simplex(tableau, cost_fn, is_maximisation);
	else solve_with_dual_simplex(tableau, cost_fn, is_maximisation);

	//If infeasible in primal, then convert problem into dual and then solve
	if (is_infeasible) {
		cout << endl << "Now converting into dual form and solving " << endl << endl;
		make_vec(cost_fn, num_vars_og);

		vvd mat1 = construct_dual_tableau(raw_matrix, cost_fn, is_maximisation, eqn_types);
		create_mat(eqn_types, mat1);
		int n_temp = mat1[0].size();
		create_maximisation_cost_function(cost_fn, is_maximisation, n_temp - 1);

		cout << "After: The coeeficient matrix: " << endl;
		printmat(mat1);
		cout << "After: The cost function: " << endl;
		pritntvecd(cost_fn);
		cout << endl;

		vvd tableau = construct_tableau(mat1);
		solve_with_dual_simplex(tableau, cost_fn, is_maximisation);
	}
}