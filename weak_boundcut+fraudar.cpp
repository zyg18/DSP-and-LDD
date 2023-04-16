#include<Eigen/Sparse>
#include<Eigen/Dense>
#include<iostream>
#include<fstream>
#include<ctime>
#include <cstdlib>
using namespace Eigen;
using namespace std;


class Mintree
{
public:
	//Mintree(int a,int b,int c,double d,int e):left(a),right(b),father(c),min_deg(d),min_source(e){}
	Mintree() { left = -1; right = -1; father = -1; min_deg = -1; min_source = -1; }
	int left;                 //Record the sequence number of the left child node in Mintree, where -1 indicates none
	int right;                //Record the serial number of the right child node in Mintree, where -1 indicates none
	int father;              //Record the serial number of the parent node in Mintree, where -1 indicates none
	double min_deg;           //Record the current maximum degree of retention
	int min_source;           //Record the source of the maximum degree currently retained, with 0 indicating left, 1 indicating right, and -1 indicating that this is the root node.
};

vector<int> kcore(SparseMatrix<double, RowMajor>& A, double bound) {
	int length = A.rows();
	vector<int>A_cut;
	int* cut_index = new int[length]();                   
	double* degree = new double[length]();
	for (int i = 0; i < length; i++) {
		degree[i] = A.row(i).sum();   
		if (degree[i] <= bound) {
			A_cut.push_back(i);
			cut_index[i] = 1;
		}
	}
	while (A_cut.size() != 0) {
		vector<int>next_cut;
		for (int i = 0; i < A_cut.size(); i++) {
			for (int j = A.outerIndexPtr()[A_cut[i]]; j < A.outerIndexPtr()[A_cut[i] + 1]; j++) {
				int num = A.innerIndexPtr()[j];
				degree[num] -= A.valuePtr()[j];
				if (degree[num] <= bound && cut_index[num] == 0) {
					next_cut.push_back(num);
					cut_index[num] = 1;
				}
			}
		}
		A_cut = next_cut;
	}
	vector<int>A_index;
	for (int i = 0; i < length; i++)
		if (cut_index[i] == 0)
			A_index.push_back(i);
	delete[]cut_index;
	delete[]degree;
	return A_index;
}


Map<SparseMatrix<double, RowMajor>>spmatrix4_weighted(SparseMatrix<double, RowMajor>& A, vector<int>A_index) {
	int length1 = A.rows();
	int length2 = A_index.size();
	int* sp_A = new int[length1]();
	for (int i = 0; i < length2; i++)
		sp_A[A_index[i]] = i + 1;                          //All indexed positions in sp_A have the subscripts of the submatrix written (to distinguish 0, the subscripts start from 1), and the remaining subscripts are 0
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


double fraudar(SparseMatrix<double, RowMajor>& A)
{
	int max_size = A.rows();
	Mintree* node = new Mintree[2 * max_size + 30];
	for (int i = 0; i < max_size; i++) {
		node[i].min_deg = A.row(i).sum();
	}
	//construct min tree
	int start_f = max_size;                                           //The starting point of the sequence number of the parent node
	int start_s = 0;                                                  //Starting point of sequence number for child nodes
	int son_size = max_size;
	while (son_size > 1) {
		int father_size = (son_size + 1) / 2;
		for (int i = 0; i < father_size; i++) {
			int left_id = start_s + 2 * i;
			int right_id = start_s + 2 * i + 1;
			int father_id = start_f + i;
			node[father_id].left = left_id;
			node[left_id].father = father_id;
			if (right_id < start_f) {
				node[father_id].right = right_id;
				node[right_id].father = father_id;
			}
			else {
				node[father_id].right = left_id;
				right_id = left_id; 
			}
			if (node[left_id].min_deg <= node[right_id].min_deg) {
				node[father_id].min_deg = node[left_id].min_deg;
				node[father_id].min_source = 0;
			}
			else {
				node[father_id].min_deg = node[right_id].min_deg;
				node[father_id].min_source = 1;
			}
		}
		son_size = father_size;
		start_s = start_f;
		start_f += father_size;
	}
	int root_id = start_s;
	double upbound = A.sum() / 2;
	double edge_sum = upbound;
	int node_num = max_size;
	double max_density = edge_sum / node_num;
	double best_edge;
	double best_node=A.rows();
	int* ensure = new int[max_size]();
	//charikar greedy algorithm
	while (node_num > 1) {                     
		edge_sum -= node[root_id].min_deg;
		node_num -= 1;
		if (edge_sum / node_num > max_density) {
			best_edge = edge_sum;
			best_node = node_num;
			max_density = edge_sum / node_num;
		}
		int del_id = (node[root_id].min_source == 0) ? node[root_id].left : node[root_id].right;
		while (node[del_id].min_source != -1)
			del_id = (node[del_id].min_source == 0) ? node[del_id].left : node[del_id].right;
		if (node[del_id].min_deg != node[root_id].min_deg) {
			cout << "when delete id:" << del_id << endl;
			cout << "node number is:" << node_num << endl;
		}
		node[del_id].min_deg = upbound;
		if (del_id < max_size && ensure[del_id] == 0)
			ensure[del_id] = 1;
		else
			cout << "something goes wrong!del_id is:" << del_id << endl;
		int father_id = node[del_id].father;
		while (father_id != -1) {
			int left_id = node[father_id].left;
			int right_id = node[father_id].right;
			if (node[left_id].min_deg <= node[right_id].min_deg) {
				node[father_id].min_source = 0;
				node[father_id].min_deg = node[left_id].min_deg;
			}
			else {
				node[father_id].min_source = 1;
				node[father_id].min_deg = node[right_id].min_deg;
			}
			father_id = node[father_id].father;
		}
		double node_degree = 0;
		for (int i = A.outerIndexPtr()[del_id]; i < A.outerIndexPtr()[del_id + 1]; i++) {
			int k = A.innerIndexPtr()[i];
			if (ensure[k] == 1)
				continue;
			else {
				double wt = A.valuePtr()[i];
				node[k].min_deg -= wt;
				int father_id = node[k].father;
				while (father_id != -1) {
					int left_id = node[father_id].left;
					int right_id = node[father_id].right;
					if (node[left_id].min_deg <= node[right_id].min_deg) {
						node[father_id].min_source = 0;
						node[father_id].min_deg = node[left_id].min_deg;
					}
					else {
						node[father_id].min_source = 1;
						node[father_id].min_deg = node[right_id].min_deg;
					}
					father_id = node[father_id].father;
				}
			}
		}
	}
	cout << "max_density is:" << max_density << "  node number is:" << best_node << endl;
	return max_density;
}




int main(int argc, char** argv)
{
	string data_in = argv[1];
	int weighted_type = atoi(argv[2]);               //0 Represents an unweighted graph, 1 represents a weighted graph, 2 represents processing as an unweighted graph, and 3 represents processing as a weighted graph
	int graph_type = atoi(argv[3]);                           //0:Monopartite graph,1：Bipartite graph
	ifstream in("./data/" + data_in);
	cout << "the data name is:" << data_in << endl;
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
				triplets.emplace_back(b - 1, a - 1, 1);                       //For symmetry
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
	A.setFromTriplets(triplets.begin(), triplets.end());         //Note that in the edgelist (4, 35) there are edges, but the sparse matrix takes 0 as the subscript, so A (3, 34) has a value of 1
	if (weighted_type == 0 || weighted_type == 2) {
		for (int i = 0; i < A.nonZeros(); i++)
			A.valuePtr()[i] = 1;
	}
	A.makeCompressed();
  double edge_sum = A.sum();
	int node_num = A.rows();
	double bound = edge_sum / (2 * node_num);
	int length = A.rows();
	vector<int>A_cut;
	int* cut_index = new int[length]();                    //Set the area that needs to be deleted to 1
	double* degree = new double[length]();
	int*ne = new int[length+1]();
	int* pr = new int[length + 1]();
	int* node_next = ne + 1;                     //node_next ：subscript range for next is -1 to length -1, node_ prev： subscript range of prev is -1 to length -1
	                                             //node_ Next [-1] indicates the first point, node_ Prev [-1] is meaningless
	int* node_prev = pr+1;
	node_next[-1] = 0;
	for (int i = 0; i < length-1; i++) {
		node_next[i] = i + 1;
		node_prev[i] = i - 1;
	}
	node_next[length - 1] = -1;
	node_prev[length - 1] = length - 2;
	clock_t start, end;
	start = clock();
	for (int i = 0; i < length; i++) {
		degree[i] = A.row(i).sum();   
		if (degree[i] < bound) {
			node_prev[node_next[i]] = node_prev[i];
			node_next[node_prev[i]] = node_next[i];
			A_cut.push_back(i);
			cut_index[i] = 1;
		}
	}
	while (A_cut.size() != 0) {
		double* degree_decrease = new double[length]();
		double edge_reduction = 0;
		for (int i = 0; i < A_cut.size(); i++) {
			edge_reduction += degree[A_cut[i]];
			for (int j = A.outerIndexPtr()[A_cut[i]]; j < A.outerIndexPtr()[A_cut[i] + 1]; j++) {
				int index = A.innerIndexPtr()[j];
				if (cut_index[index] == 0)
					degree_decrease[index] += A.valuePtr()[j];
			}
		}
		for (int i = node_next[-1]; i!=-1; i=node_next[i]){
			degree[i] -= degree_decrease[i];
			edge_reduction += degree_decrease[i];
		}
		node_num -= A_cut.size();
		edge_sum -= edge_reduction;
		bound = edge_sum / (2 * node_num);
		A_cut.clear();
		for (int i = node_next[-1]; i != -1; i = node_next[i]) {
			if (degree[i] < bound) {
				A_cut.push_back(i);
				cut_index[i] = 1;
				node_prev[node_next[i]] = node_prev[i];
				node_next[node_prev[i]] = node_next[i];
			}
		}
	}
	vector<int>A_index;
	for (int i = 0; i < length; i++)
		if (cut_index[i] == 0)
			A_index.push_back(i);
	delete[]cut_index;
	delete[]degree;
	A = spmatrix4_weighted(A, A_index);
 
  
	fraudar(A);
	end = clock();
	double time = (double)(end - start) / CLOCKS_PER_SEC;
	cout << "weak_boundcut+fraudar cost time is:" << time << "s" << endl;
	return 0;
}