#include<Eigen/Sparse>
#include<Eigen/Dense>
#include<iostream>
#include<fstream>
#include<ctime>
#include<iomanip>
#include <cstdlib>
using namespace Eigen;
using namespace std;



Map<SparseMatrix<double, RowMajor>>spmatrix4_weighted(SparseMatrix<double, RowMajor>& A, vector<int>A_index) {
	int length1 = A.rows();
	int length2 = A_index.size();
	int* sp_A = new int[length1]();
	for (int i = 0; i < length2; i++)
		sp_A[A_index[i]] = i + 1;                           
	int* indptr = new int[length2 + 1];
	indptr[0] = 0;
	vector<int>indice;
	vector<double>data;
	int ptr = 0;
	for (int i = 0; i < length2; i++) {
		int start = A.outerIndexPtr()[A_index[i]];
		int fin = A.outerIndexPtr()[A_index[i] + 1];
		for (int k = start; k < fin; k++) {
			int ind = sp_A[A.innerIndexPtr()[k]];
			if (ind != 0) {
				indice.push_back(ind - 1);
				data.push_back(A.valuePtr()[k]);
				ptr += 1;
			}
		}
		indptr[i + 1] = ptr;
	}
	int* arr1 = new int[indice.size()];
	std::copy(indice.begin(), indice.end(), arr1);
	double* arr2 = new double[data.size()];
	std::copy(data.begin(), data.end(), arr2);
	Map<SparseMatrix<double, RowMajor>>sm1(length2, length2, ptr, indptr, arr1, arr2);
	return sm1;
}



struct Edge {
	int y, next;
};

struct Node {                                
	int deg, next, prev, idx;
	inline void clear() {
		deg = next = prev = 0;
		idx = -1;
	}
};

Node* lists;

__inline void linklists(int x, int y) {     
	if (y == 0) return;
	lists[x].next = y;
	lists[y].prev = x;
};

int* nxt, * prv, * itr;                      

__inline void linknodes(int x, int y) {     
	if (y == -1) return;
	nxt[x] = y;
	prv[y] = x;
};

__inline void eraselist(int x) {            
	lists[lists[x].prev].next = lists[x].next;
	if (lists[x].next != 0) lists[lists[x].next].prev = lists[x].prev;
};

__inline void erasenode(int x) {           
	if (prv[x] == -1) {
		lists[itr[x]].idx = nxt[x];
	}
	if (prv[x] != -1) nxt[prv[x]] = nxt[x];
	if (nxt[x] != -1) prv[nxt[x]] = prv[x];

};

int l = 0;
Edge* edges;
int* idx;

__inline void build(int x, int y) {         
	edges[++l].next = idx[x];
	edges[l].y = y;
	idx[x] = l;
};


int main(int argc, char** argv) {
	string data_in=argv[1];
  int weighted_type=atoi(argv[2]);               
  int graph_type=atoi(argv[3]);                           
	ifstream in("./data/"+data_in);
  cout<<endl<<"the data name is:"<<data_in<<endl;
	if (!in) {
		cout << "fail to open the file" << endl;
		return -1;
	}
	else
		cout << "opened the file" << endl;
	string comment;
	getline(in, comment);									
	char s;
	in >> s;
	int edge_num, m, n;									    
	in >> edge_num >> m >> n;
	cout << "edge_num=" << edge_num << "  " << "m=" << m << "  " << "n=" << n << "  " << endl;
	int max_m = 0, max_n = 0;
	int max_size = 0;
	vector<Triplet<int>>triplets;
	if (graph_type == 0) {
		if (weighted_type == 0) {
			int a, b;
			while (in >> a >> b)                  
			{
				if (a == b)
					continue;
				triplets.emplace_back(a - 1, b - 1, 1);
				triplets.emplace_back(b - 1, a - 1, 1);                       
				max_m = max_m > a ? max_m : a;
				max_n = max_n > b ? max_n : b;
			}
		}
		else if (weighted_type == 1) {
			int a, b, c;
			while (in >> a >> b >> c)                  
			{
				if (a == b)
					continue;
				triplets.emplace_back(a - 1, b - 1, c);
				triplets.emplace_back(b - 1, a - 1, c);                       
				max_m = max_m > a ? max_m : a;
				max_n = max_n > b ? max_n : b;
			}
		}
		else if (weighted_type == 2) {          
			int a, b, c;
			while (in >> a >> b >> c)
			{
				if (a == b)
					continue;
				triplets.emplace_back(a - 1, b - 1, 1);
				triplets.emplace_back(b - 1, a - 1, 1);                       
				max_m = max_m > a ? max_m : a;
				max_n = max_n > b ? max_n : b;
			}
		}
		else if (weighted_type == 3) {          
			int a, b, c, d;
			while (in >> a >> b >> c >> d)
			{
				if (a == b)
					continue;
				triplets.emplace_back(a - 1, b - 1, c);
				triplets.emplace_back(b - 1, a - 1, c);                       
				max_m = max_m > a ? max_m : a;
				max_n = max_n > b ? max_n : b;
			}
		}
		else {
			cout << "Unrecognized type" << endl;
			return 0;
		}
		max_size = max_m > max_n ? max_m : max_n;
	}
	else if (graph_type == 1) {
		if (weighted_type == 0) {
			int a, b;
			while (in >> a >> b)                                           
			{
				max_m = max_m > a ? max_m : a;
				max_n = max_n > b ? max_n : b;
				b = m + b;
				triplets.emplace_back(a - 1, b - 1, 1);
				triplets.emplace_back(b - 1, a - 1, 1);                       
			}
		}
		else if (weighted_type == 1) {
			int a, b, c;
			while (in >> a >> b >> c)                 
			{
				max_m = max_m > a ? max_m : a;
				max_n = max_n > b ? max_n : b;
				b = m + b;
				triplets.emplace_back(a - 1, b - 1, c);
				triplets.emplace_back(b - 1, a - 1, c);                       

			}
		}
		else if (weighted_type == 2) {          
			int a, b, c;
			while (in >> a >> b >> c)
			{
				max_m = max_m > a ? max_m : a;
				max_n = max_n > b ? max_n : b;
				b = m + b;
				triplets.emplace_back(a - 1, b - 1, 1);
				triplets.emplace_back(b - 1, a - 1, 1);                      
			}
		}
		else if (weighted_type == 3) {          
			int a, b, c, d;
			while (in >> a >> b >> c >> d)
			{
				max_m = max_m > a ? max_m : a;
				max_n = max_n > b ? max_n : b;
				b = m + b;
				triplets.emplace_back(a - 1, b - 1, c);
				triplets.emplace_back(b - 1, a - 1, c);                      
			}
		}
		else {
			cout << "Unrecognized type" << endl;
			return 0;
		}
		max_size = m + n;
	}
	cout << "max_m is:" << max_m << "  " << "max_n is:" << max_n << endl;
	SparseMatrix<double, RowMajor>A(max_size, max_size);                       
	A.setFromTriplets(triplets.begin(), triplets.end());         
	if (weighted_type == 0 || weighted_type == 2) {
		for (int i = 0; i < A.nonZeros(); i++)
			A.valuePtr()[i] = 1;
	}
	A.makeCompressed();
	int edge_sum = A.sum();
	int node_num = A.rows();
	clock_t start, end;
	start = clock();
	double bound = (double)edge_sum / (2 * node_num);
	int length = A.rows();
	vector<int>A_cut;
	int* cut_index = new int[length]();                    
	int* degree = new int[length]();
	int max_degree=0;
	for (int i = 0; i < length; i++) {
		degree[i] = A.row(i).sum();   //第i个点的度
		max_degree = max_degree > degree[i] ? max_degree:degree[i];
		if (degree[i] < bound) {
			A_cut.push_back(i);
			cut_index[i] = 1;
		}
	}
	Node* node = new Node[length];
	int* degreelist = new int[max_degree+1]();
	for (int i = 0; i <= max_degree; i++)
		degreelist[i] = -1;
	//建立节点列表，列表内部是双向链表
	for (int i = 0; i <length; i++) {
		node[i].next = degreelist[degree[i]];
		node[i].prev = -1;
		if (degreelist[degree[i]] != -1)
			node[degreelist[degree[i]]].prev= i;
		degreelist[degree[i]] = i;
	}
	int iters = 1;
	int*is_dealed = new int[length]();                         
	while (A_cut.size() != 0) {
		int* degree_decrease = new int[length]();
		int edge_reduction = 0;
		vector<int>neighbour;
		for (int i = 0; i < A_cut.size(); i++) {
			edge_reduction += degree[A_cut[i]];
			for (int j = A.outerIndexPtr()[A_cut[i]]; j < A.outerIndexPtr()[A_cut[i] + 1]; j++) {
				int index = A.innerIndexPtr()[j];
				if (cut_index[index] == 0) {
					degree_decrease[index] += 1;
					if (is_dealed[index] == iters)
						continue;
					is_dealed[index] = iters;
					neighbour.push_back(index);
				}
			}
		}
		for (int i = 0; i < neighbour.size(); i++) {
			degree[neighbour[i]] -= degree_decrease[neighbour[i]];
			edge_reduction += degree_decrease[neighbour[i]];
		}
		int deg = floor(bound);
		degreelist[deg] = -1;
		node_num -= A_cut.size();
		edge_sum -= edge_reduction;
		bound = (double)edge_sum / (2 * node_num);
		A_cut.clear();
		for (int i = 0; i < neighbour.size(); i++) {
			int index = neighbour[i];
			if (node[index].prev != -1)
				node[node[index].prev].next = node[index].next;
			else
				degreelist[degree[index] + degree_decrease[index]] = node[index].next;
			if (node[index].next != -1)
				node[node[index].next].prev = node[index].prev;
			if (degree[index] > bound) {
				node[index].next = degreelist[degree[index]];
				node[index].prev = -1;
				if (degreelist[degree[index]] != -1)
					node[degreelist[degree[index]]].prev = index;
				degreelist[degree[index]] = index;
			}
			else {
				node[index].next = degreelist[deg];
				node[index].prev = -1;
				if (degreelist[deg] != -1)
					node[degreelist[deg]].prev = index;
				degreelist[deg] = index;
			}
		}
		for (int i = deg; i <= floor(bound); i++) {
			int index = degreelist[i];
			while (index != -1) {
				cut_index[index] = 1;
				A_cut.push_back(index);
				index = node[index].next;
			}
		}
		iters++;
	}
	cout << "bound is:" << bound << " node_num is:" << node_num << endl;
	end = clock();
  vector<int>A_index;
	for (int i = 0; i < length; i++)
		if (cut_index[i] == 0)
			A_index.push_back(i);
	delete[]cut_index;
	delete[]degree;
	double time = (double)(end - start) / CLOCKS_PER_SEC;
	cout << "weak_boundcut cost time is:" << time << "s" << endl;

	//unweighted greedy++ initialization
  A = spmatrix4_weighted(A, A_index);
	int node_number = A.rows();
	int edge_number = A.nonZeros()/2;
	edges = new Edge[edge_number * 2 + 10];
	idx = new int[node_number];
	memset(idx, 0, sizeof(int)* node_number);
	int* init_deg = new int[node_number];
	memset(init_deg, 0, sizeof(int)* node_number);
	l = 0;
	lists = new Node[node_number + 2 * edge_number + 10];
	int n_list = 0;
	itr = new int[node_number];
	int* deg = new int[node_number], * w = new int[node_number], * pos = new int[node_number];
	memset(deg, 0, sizeof(int)* node_number);
	memset(w, 0, sizeof(int)* node_number); //initial vertex weights=0, i.e., no self loops at the start 
	memset(pos, 0, sizeof(int)* node_number);
	prv = new int[node_number]; nxt = new int[node_number];
	for (int i = 0; i < node_number; i++) {
		int start = A.outerIndexPtr()[i];
		int fin = A.outerIndexPtr()[i + 1];
		for (int j = start; j < fin; j++) {
			int k = A.innerIndexPtr()[j];
			if (i >= k)
				continue;
			build(i, k);
			build(k, i);
			init_deg[i]++;
			init_deg[k]++;
		}
	}
	pair<int, int>* deg_sorted = new pair<int, int>[node_number];
	vector<int> m_ans;
	double mm_density = 0;
	clock_t end2 = clock();
	time = (double)(end2 - end) / CLOCKS_PER_SEC;
	cout << "unweighted greedy++ initialization cost time is:" << time << "s" << endl;



  ofstream output;
  output.open("./output/uwgreedy_"+data_in+".csv",ios::out|ios::trunc); 
	clock_t run_begin, run_end;
	iters = 2000;
	double sum_iter_times = 0;
	for (int tt = 0; tt < iters; tt++) {
		run_begin = clock();
		for (int i = 0; i < node_number; i++) {
			nxt[i] = prv[i] = -1;
			pos[i] = 0;
			deg[i] = w[i] + init_deg[i];                                       //degree for this iteration is "vertex weight" + actual degree
			deg_sorted[i] = make_pair(deg[i], i);
		}
		sort(deg_sorted, deg_sorted + node_number);                                    
		n_list = 0;
		for (int i = 0; i < node_number; i++) {
			int v = deg_sorted[i].second;                                     
			if (n_list == 0 || lists[n_list].deg != deg_sorted[i].first) {     
				++n_list;
				lists[n_list].clear();
				linklists(n_list - 1, n_list);
				lists[n_list].deg = deg_sorted[i].first;
			}
			linknodes(v, lists[n_list].idx);
			lists[n_list].idx = v;
			itr[v] = n_list;
		}
		int edge_sum = A.nonZeros()/2;
		double max_density = (double)edge_sum / node_number;
		int cur_m = edge_sum, cur_n = node_number;
		//vector<int> ans;
		int max_size = A.rows();
		while (lists[0].next) {
			int i = lists[0].next;
			int k = lists[i].idx;
			if (nxt[k] == -1) {
				eraselist(i);
			}
			else {
				erasenode(k);
			}
			pos[k] = -1;
			w[k] = deg[k]; //increment vertex weight for the next iteration (self loops)
			cur_n -= 1;
			//ans.push_back(k);
			for (int p = idx[k]; p; p = edges[p].next) { //decrement degrees of k's neighbors
				int j = edges[p].y;
				if (pos[j] == -1) continue;
				cur_m -= 1;
				int i = itr[j];
				erasenode(j);
				int i1 = lists[i].prev;
				if (lists[i].idx == -1) eraselist(i);
				deg[j]--;
				prv[j] = nxt[j] = -1;
				if (i1 == 0 || lists[i1].deg != deg[j]) {
					++n_list;
					lists[n_list].clear();
					itr[j] = n_list;
					int i2 = lists[i1].next;
					lists[n_list].deg = deg[j];
					lists[n_list].idx = j;
					linklists(i1, n_list);
					if (i2) linklists(n_list, i2);
				}
				else {
					linknodes(j, lists[i1].idx);
					lists[i1].idx = j;
					itr[j] = i1;
				}
			}
			if (cur_n == 0) continue;
			if (max_density < (double)cur_m / cur_n) {
				max_size = cur_n;
			}
			max_density = max(max_density, (double)cur_m / cur_n);
		}
		//reverse(ans.begin(), ans.end());
		//ans.resize(node_number - max_size);
		if (max_density > mm_density) {
			//m_ans = ans;
			mm_density = max_density;
		}
		double run_time;
		run_end = clock();
		run_time = (double)(run_end - run_begin) / CLOCKS_PER_SEC;
		sum_iter_times += run_time;
		cout << "Max size is:" << max_size << endl;
		cout << "Max density until iteration " << tt + 1 << ": " << mm_density << endl;
		cout << "Avg time per iteration: " << sum_iter_times / (tt + 1) << "s" << endl;
    output<<tt<<","<<setprecision(8)<<mm_density<<","<<sum_iter_times<<","<<max_size<<endl;
	}
	//end = clock();
	//time = (double)(end - start) / CLOCKS_PER_SEC;
	//cout << "boundcut cost time is:" << time << "s" << endl;
	return 0;
}
