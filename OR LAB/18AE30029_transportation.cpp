//Lab09 revided submission
//Transportation problem

//Note: The inputs for the questions and an example have been given at the end of the program in the form of comment for your convenience

#include<bits/stdc++.h>
using namespace std;
#define vi vector<int>
#define vvi vector<vector<int>>
#define pi pair<int,int>
#define vb vector<bool>
#define vvb vector<vector<bool> >
#define vpi vector<pair<int,int>>
#define num_iterations_max 10

int num_origins;
int num_destinations;

//Pretty print table
void print_mat(vvi mat){
	for(int i=0;i<mat.size();i++){
		for(int j=0;j<mat[i].size();j++){
			cout.width(5);
			cout<<mat[i][j]<<" ";
		}
		cout<<endl;
	}
	cout<<endl;
}

void printvec(vi arr){
	for(int i=0;i<arr.size();i++)
		cout<<arr[i]<<" ";
	cout<<endl;
}

vvi construct_talbeau(){
	cout<<"Enter the number of origins: ";
	cin>>num_origins;
	cout<<"Enter the number of destinations: ";
	cin>>num_destinations;

	vvi tableau(num_origins+1 , vector<int> (num_destinations+1, 0));

	cout<<endl<<"Enter the supplies at the "<<num_origins<<" origins:"<<endl;
	for(int i=0;i<num_origins;i++){
		cin>>tableau[i][num_destinations];
	}

	cout<<endl<<"Enter the demands at the "<<num_destinations<<" destinations:"<<endl;
	for(int i=0;i<num_destinations;i++){
		cin>>tableau[num_origins][i];
	}

	cout<<endl<<"Enter the transportation costs from source i to destination j in row wise order"<<endl;
	for(int i=0;i<num_origins;i++){
		for(int j=0;j<num_destinations;j++){
			cin>>tableau[i][j];
		}
	}

	int total_supply = 0;
	for(int i=0;i<num_origins;i++){
		total_supply+=tableau[i][num_destinations];
	}
	tableau[num_origins][num_destinations] = total_supply;

	cout<<endl<<"The initial problem: "<<endl;
	return tableau;
}

bool is_TP_balanced(vvi tableau){
	int total_supply = 0;
	int total_demand = 0;
	for(int i=0;i<num_origins;i++){
		total_supply+=tableau[i][num_destinations];
	}
	for(int i=0;i<num_destinations;i++){
		total_demand+=tableau[num_origins][i];
	}
	if(total_demand == total_supply) return 1;
	return 0;
}

vvi NW_corner_method(vvi tableau){
	int curr_row = 0;
	int curr_col = 0;
	vvi mat(num_origins , vector<int> (num_destinations, -1));

	int curr_demand = tableau[num_origins][0];
	int curr_supply = tableau[0][num_destinations];

	while((curr_row < num_origins) && (curr_col < num_destinations)){
		if(curr_supply <= curr_demand){//go down
			mat[curr_row++][curr_col] = curr_supply;
			curr_demand-=curr_supply;
			if(curr_row >= num_origins) break;
			curr_supply=tableau[curr_row][num_destinations];
		}
		else{//go right
			mat[curr_row][curr_col++] = curr_demand;
			curr_supply-=curr_demand;
			if(curr_col >= num_destinations) break;
			curr_demand=tableau[num_origins][curr_col];
		}
	}
	return mat;
}

int find_first_null_ui(vb Ui_bool_vec){
	for(int i=0;i<Ui_bool_vec.size();i++){
		if(Ui_bool_vec[i]==0) return i;
	}
	return -1;
}


// Ui+Vj = Cij ... calculate Ui, Vj from the allocated cells
void calculate_UVs(vvi &mat, vvi tableau, vi &Ui_vec, vi &Vj_vec){
	//cout<<"/////CALCULATING UIS"<<endl;
	int i, num_Uis_found, num_Vjs_found;
	i = num_Vjs_found = num_Uis_found = 0;

	vb Ui_bool_vec(Ui_vec.size(), 0);
	vb Vj_bool_vec(Vj_vec.size(), 0);

	//start with U1=0 always
	Ui_vec[0]=0; num_Uis_found++;
	Ui_bool_vec[0] = 1;

	while(1){
		int num_iters_max = 20;
		int iter_index = 0;
		bool flag_break=0;
		while((num_Uis_found < num_origins) || (num_Vjs_found < num_destinations)){
			iter_index++;
			if(iter_index > num_iters_max){
				flag_break = 1;
				break;
			}
			for(int i=0;i<num_origins;i++){
				for(int k=0;k<2;k++){
					for(int j=0;j<num_destinations;j++){
						//cout<<"i,j: "<<i<<","<<j<<endl;
						if(mat[i][j]!=-1){//Only consider allocated cells
							if((Ui_bool_vec[i]) && (Vj_bool_vec[j])) continue;
							else if(Ui_bool_vec[i]){
								Vj_vec[j] = tableau[i][j] - Ui_vec[i];//update Vj_vec[j]
								//cout<<"vj["<<j<<"] updated to "<<Vj_vec[j]<<endl;
								Vj_bool_vec[j]=1;
								num_Vjs_found++;
							}
							else if(Vj_bool_vec[j]){
								Ui_vec[i] = tableau[i][j] - Vj_vec[j];//update Ui_vec[i]
								//cout<<"ui["<<i<<"] updated to "<<Ui_vec[i]<<endl;
								Ui_bool_vec[i]=1;
								num_Uis_found++;
							}
						}
					}
				}
			}
		}

		if(flag_break){
			int i = find_first_null_ui(Ui_bool_vec);
			if(i==-1){
				cout<<"Kuch laphda hai, how are all Uis filled and still infiinte loop!?"<<endl;
			}
			Ui_bool_vec[i]=1;
			Ui_vec[i]=0;
			num_Uis_found++;
		}
		else break;
	}

	//cout<<"/////DONE CALCULATNIG UIS"<<endl;
}

// Pij = Ui + Vj - Cij ... calculate penalties for unallocated cells
// If all Pij are <=0 means terminate
// Else choose (i,j) from the greatest Pij as the new basic cell(NBC)
void calculate_NBC(vvi &mat, vvi tableau, vi Ui_vec, vi Vj_vec, pi &nbc){
	cout<<"Penalties: "<<endl;
	int max_till_now = -1;
	for(int i=0;i<num_origins;i++){
		for(int j=0;j<num_destinations;j++){
			if(mat[i][j]==-1){//Only  consider Unallocated cells
				int curr_Pij = Ui_vec[i] + Vj_vec[j] - tableau[i][j];
				cout<<"P["<<i<<"]"<<"["<<j<<"]: "<<curr_Pij<<endl;
				if(curr_Pij > 0){
					if(curr_Pij > max_till_now){
						max_till_now = curr_Pij;
						nbc.first = i; nbc.second = j;
					}
				}
			}
		}
	}
}

bool is_valid_cell(pi cell){
	if((cell.first < num_origins) && (cell.second < num_destinations) && (cell.first >= 0) && (cell.second >= 0)) return 1;
	return 0;
}

bool are_in_same_line(pi prev, pi curr, pi next){
	if((prev.first == curr.first) && (next.first == curr.first)) return 1;
	if((prev.second == curr.second) && (next.second == curr.second)) return 1;
	return 0;
}

//Modify mat to be used in the next iteration to calculate new values of Ui, Uj
//Make a closed loop consisting of only horizontal, vertical lines such that the loop
//Passes only through allocated cells
//Use DFS(depth first search) to implement the LOOP detection
bool is_loop_found_func(vvi mat, vpi next_cell_options, pi nbc, pi prev_cell, pi curr_cell){
	for(int i=0;i<next_cell_options.size();i++){
		if((next_cell_options[i].first == nbc.first) && (next_cell_options[i].second == nbc.second)){
			if((mat[curr_cell.first][curr_cell.second]!=-1) || (are_in_same_line(prev_cell, curr_cell, next_cell_options[i])))
				return 1;
		}
	}
	return 0;
}

void printstack(stack<pi> st){
	cout<<"STACK-> ";
	while(!st.empty()){
		pi top = st.top();
		cout<<"("<<top.first<<","<<top.second<<"), ";
		st.pop();
	}
	cout<<endl;
}

void find_loop(stack<pi> &loop, vvi mat, vvi tableau, pi prev_cell, pi curr_cell, vvb &flag_mat, bool &is_loop_found, pi nbc, int &length_of_cycle){
	vpi next_cell_options = {{curr_cell.first+1, curr_cell.second}, {curr_cell.first, curr_cell.second+1}, {curr_cell.first-1, curr_cell.second}, {curr_cell.first, curr_cell.second-1}};
	if(is_loop_found) return;
	if(is_loop_found_func(mat, next_cell_options, nbc, prev_cell, curr_cell)){
		if((length_of_cycle > 2)){
			cout<<"Found a loop"<<endl;
			is_loop_found = 1;
			return;
		}
	}

	for(int i=0;i<next_cell_options.size();i++){
		pi next_cell = next_cell_options[i];
		//cout<<"cell: ("<<next_cell.first<<","<<next_cell.second<<"): ";
		if(!is_valid_cell(next_cell)){
			//cout<<"NOT VALID"<<endl;
			continue;
		}
		if(flag_mat[next_cell.first][next_cell.second] == 0){
			if(are_in_same_line(prev_cell, curr_cell, next_cell) || (mat[curr_cell.first][curr_cell.second]!=-1)){
				loop.push(next_cell);
				//cout<<"("<<next_cell.first<<","<<next_cell.second<<") pushed, loc: "<<length_of_cycle<<", "; printstack(loop);
				flag_mat[next_cell.first][next_cell.second] = 1;
				length_of_cycle++;
				find_loop(loop, mat, tableau, curr_cell, next_cell, flag_mat, is_loop_found, nbc, length_of_cycle);
				if(is_loop_found) return;
				flag_mat[next_cell.first][next_cell.second] = 0;
				loop.pop();
				//cout<<"("<<next_cell.first<<","<<next_cell.second<<") popped, loc: "<<length_of_cycle<<", "; printstack(loop);
				length_of_cycle--;
			}
			else{
				//cout<<"Not in same line"<<endl;
			}
		}
		else{
			//cout<<"already visited"<<endl;
		}
	}
}

vpi form_turn_index_list(vpi loop_vec){
	vpi ret_vec;
	int n = loop_vec.size();

	for(int i=0;i<n;i++){
		pi curr = loop_vec[i];
		pi next = loop_vec[(i+1)%n];
		pi prev = loop_vec[(i-1+n)%n];
		if(!are_in_same_line(prev, curr, next)){
			ret_vec.push_back(loop_vec[i]);
		}
	}
	return ret_vec;
}

void modify_mat(vvi &mat, vvi tableau, pi nbc, bool &is_loop_found){
	vvb flag_mat(num_origins , vector<bool> (num_destinations, 0));//Is the cell already visited?
	flag_mat[nbc.first][nbc.second]=1;

	stack<pi> loop;
	loop.push(nbc);
	//cout<<endl<<endl<<endl<<"nbc("<<nbc.first<<","<<nbc.second<<") pushed, "; printstack(loop);

	int length_of_cycle = 1;
	find_loop(loop, mat, tableau, nbc, nbc, flag_mat, is_loop_found, nbc, length_of_cycle);
	if(!is_loop_found){
		cout<<endl<<"Couldn't find loop,  ";
		return;
	}

	vpi loop_vec;
	while(!loop.empty()){
		pi top = loop.top();
		//cout<<"("<<top.first<<","<<top.second<<"), ";
		//loop_vec.push_back(top);
		loop_vec.insert(loop_vec.begin(), top);
		loop.pop();
	}
	cout<<"The detected loop: ";
	for(int i=0;i<loop_vec.size();i++){
		cout<<"("<<loop_vec[i].first<<","<<loop_vec[i].second<<"), ";
	}
	cout<<endl;


	vpi turn_indices = form_turn_index_list(loop_vec);
	//Now that we have detected the turn_indices, now update mat
	//Only change the elements where the loop is turning
	int min_till_now = INT_MAX;
	for(int i=0;i<turn_indices.size();i++){
		if(i%2){
			pi curr_cell = turn_indices[i];
			if(mat[curr_cell.first][curr_cell.second] == -1){
				cout<<"LAPHDA since turning point cant be unallocated"<<endl;
				return;
			}
			if(mat[curr_cell.first][curr_cell.second] < min_till_now){
				min_till_now = mat[curr_cell.first][curr_cell.second];
			}
		}
	}

	for(int i=0;i<turn_indices.size();i++){
		pi curr_cell = turn_indices[i];
		if(i%2){
			if(mat[curr_cell.first][curr_cell.second] - min_till_now == 0)
				mat[curr_cell.first][curr_cell.second] = -1;
			else
				mat[curr_cell.first][curr_cell.second] -= min_till_now;
		}
		else{
			if(mat[curr_cell.first][curr_cell.second] == -1)
				mat[curr_cell.first][curr_cell.second] = min_till_now;
			else
				mat[curr_cell.first][curr_cell.second] += min_till_now;
		}
	}
}

pi find_min_cell(vvi tableau, vb row_taken, vb col_taken){
	int min = INT_MAX;
	pi min_pair = {-1,-1};
	for(int i=0;i<num_origins;i++){
		for(int j=0;j<num_destinations;j++){
			if(!row_taken[i] && !col_taken[j]){
				if(tableau[i][j] < min){
					min_pair = {i,j};
					min = tableau[i][j];
				}
			}
		}
	}
	return min_pair;
}

//Returns mat calculated using the least cost method
vvi LeastCostMethod(vvi tableau){
	int curr_row = 0;
	int curr_col = 0;
	vvi mat(num_origins , vector<int> (num_destinations, -1));
	vb row_taken(num_origins,0);
	vb col_taken(num_destinations,0);

	while(1){
		pi cell_min = find_min_cell(tableau, row_taken, col_taken);
		if(cell_min == make_pair(-1,-1)) break;

		//cout<<"cell_min: ("<<cell_min.first<<","<<cell_min.second<<")"<<endl;

		if(tableau[cell_min.first][num_destinations] < tableau[num_origins][cell_min.second]){//supply < demand
			row_taken[cell_min.first] = 1;
			mat[cell_min.first][cell_min.second] = tableau[cell_min.first][num_destinations];
			tableau[num_origins][cell_min.second] -= tableau[cell_min.first][num_destinations];
			tableau[cell_min.first][num_destinations]=0;
		}
		else{//demand < supply
			col_taken[cell_min.second] = 1;
			mat[cell_min.first][cell_min.second] = tableau[num_origins][cell_min.second];
			tableau[cell_min.first][num_destinations] -= tableau[num_origins][cell_min.second];
			tableau[num_origins][cell_min.second]=0;
		}
	}

	return mat;
}

int find_cost(vvi mat, vvi tableau){
	int cost=0;
	for(int i=0;i<mat.size();i++){
		for(int j=0;j<mat[i].size();j++){
			if(mat[i][j]!=-1){
				cost += (mat[i][j]*tableau[i][j]);
			}
		}
	}
	return cost;
}

void solve_using_uv_method(vvi tableau, vvi mat, bool method){
	int iter_index = 1;
	while(1){
		cout<<endl<<"Iteration: "<<iter_index<<endl;
		cout<<"Current cost: "<<find_cost(mat, tableau)<<endl;
		vi Ui_vec(num_origins, -1);
		vi Vj_vec(num_destinations, -1);

		calculate_UVs(mat, tableau, Ui_vec, Vj_vec);
		cout<<"Ui_vec: "; printvec(Ui_vec);
		cout<<"Vj_vec: "; printvec(Vj_vec);

		pi nbc{-1,-1};//New Basic Cell -> Will remain (-1,-1) if optimality has been reached
		calculate_NBC(mat, tableau, Ui_vec, Vj_vec, nbc);
		if((nbc.first == -1) && (nbc.second == -1)){
			cout<<"optimality has been reached!"<<endl;
			break;
		}
		cout<<"New Basic Cell: ("<<nbc.first<<", "<<nbc.second<<")"<<endl;

		bool is_loop_found = 0;
		modify_mat(mat, tableau, nbc, is_loop_found);
		if(!is_loop_found){
			if(!method){
				cout<<"Moving on to try The Least Cost Method"<<endl;
				cout<<endl<<endl<<"Least Cost Method: "<<endl;
				mat = LeastCostMethod(tableau);
				cout<<endl<<"Initial table(-1 represent unallocated cells) -> "<<endl;print_mat(mat);

				method=1;//LeastCostMethod
				solve_using_uv_method(tableau, mat, method);
				return;
			}
			else{
				cout<<"Both initialisation methods failed....RETURNING!!"<<endl;
				return;
			}
		}

		cout<<"Table "<<iter_index<<" -> "<<endl;
		print_mat(mat);
		//print_mat(tableau);
		if(iter_index > num_iterations_max){
			cout<<"Iterations exceeded... LAPHDA.... RETURNING"<<endl;
			break;
		}
		iter_index++;
	}
	cout<<endl<<endl<<"Final cost: "<<find_cost(mat, tableau)<<endl;
}

vvi balance_tableau(vvi tableau){
	int supply = 0;
	int demand = 0;
	for(int i=0;i<num_origins;i++){
		supply += tableau[i][num_destinations];
	}
	for(int j=0;j<num_destinations;j++){
		demand += tableau[num_origins][j];
	}
	cout<<"Total supply: "<<supply<<", Total demand: "<<demand<<endl;
	int num_rows, num_cols;
	if(demand < supply){//dummy demand center
		num_rows = num_origins+1;
		num_cols = num_destinations+2;
	}
	else if(demand > supply){//dummy supply center
		num_rows = num_origins+2;
		num_cols = num_destinations+1;
	}
	else{
		cout<<"ABEY <Balanced hi to hai> !!"<<endl;
	}


	vvi ret_tableau(num_rows , vector<int> (num_cols, 0));
	for(int i=0;i<num_origins+1;i++){
		for(int j=0;j<num_destinations+1;j++){
			ret_tableau[i][j] = tableau[i][j];
		}
	}

	if(num_rows == num_origins+2){
		for(int j=0;j<num_destinations;j++){
			ret_tableau[num_origins][j] = 0;
			ret_tableau[num_origins+1][j] = tableau[num_origins][j];
		}

	}
	else if(num_cols == num_destinations+2){
		for(int i=0;i<num_origins;i++){
			ret_tableau[i][num_destinations] = 0;
			ret_tableau[i][num_destinations+1] = tableau[i][num_destinations];
		}
	}

	num_origins = num_rows-1;
	num_destinations = num_cols-1;
	return ret_tableau;
}

int main(){
	vvi tableau = construct_talbeau();
	bool is_balanced = is_TP_balanced(tableau);
	if(!is_balanced){
		cout<<"The given transportation problem is not balanced. Adding dummy variables to balance it!"<<endl;
		//Modify tableau by adding dummy columns/rows to it

		//Whenever total demand exceeds total supply, we introduce a dummy supply centre
		//(additional supply row) in the transportation problem to meet extra demand.
		//The unit transportation cost for the cells of its dummy row are set equal to zero.
		tableau = balance_tableau(tableau);
	}
	print_mat(tableau);

	//cout<<"The given transportation problem is a Balanced TP!"<<endl;
	cout<<"Finding Initital BFS using NorthWest Corner Method"<<endl;

	vvi mat = NW_corner_method(tableau);
	cout<<endl<<"Initial table(-1 represent unallocated cells) -> "<<endl;
	print_mat(mat);
	//print_mat(tableau);

	bool method=0;//NW_corner_method
	solve_using_uv_method(tableau, mat, method);

	return 0;
}


/*
Example problem
3
4
250 350 400
200 300 350 150
3 1 7 4
2 6 5 9
8 3 3 2
Q1
3
4
3 5 7
3 2 6 4
10 7 3 6
1 6 8 3
7 4 5 3
Q2
3
4
7 9 18
5 8 7 14
19 30 50 10
70 30 40 60
40 8 70 20
Q3
3
4
30 50 80
20 60 55 40
3 8 7 4
5 2 9 5
4 3 6 2
*/
