#include<Eigen/Sparse>
#include<Eigen/Dense>
#include<iostream>
#include<fstream>
#include<ctime>
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




int main(int argc, char** argv)
{
	string data_in = argv[1];
	int weighted_type = atoi(argv[2]);               
	int graph_type = atoi(argv[3]);                           
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
		max_size = max_m + max_n;
	}
	cout << "max_m is:" << max_m << "  " << "max_n is:" << max_n << endl;
	SparseMatrix<double, RowMajor>A(max_size, max_size);                       
	A.setFromTriplets(triplets.begin(), triplets.end());         
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
	int* cut_index = new int[length]();                    
	double* degree = new double[length]();
	int*ne = new int[length+1]();
	int* pr = new int[length + 1]();
	int* node_next = ne + 1;                     
	                                             
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
    cout<<A_cut.size()<<endl;
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
  A=spmatrix4_weighted(A,A_index);
	delete[]cut_index;
	delete[]degree;
	cout << "bound is:" << bound << " node_num is:" << node_num << endl;
	end = clock();
	double time = (double)(end - start) / CLOCKS_PER_SEC;
	cout << "weak_boundcut cost time is:" << time << "s" << endl;
	cout << endl;
  ofstream p;
  p.open("./FinalCodeAndAppendix CameraReady/datasets_cleaned/pruning_"+data_in,ios::out|ios::trunc);
  length=A.rows();
  p<<length<<" "<<A.nonZeros()/2<<endl;
  for(int i=0;i<length; i++)
  		for (int j = A.outerIndexPtr()[i]; j < A.outerIndexPtr()[i + 1]; j++) {
        int index = A.innerIndexPtr()[j];
        if(index<i)
          continue;
        p<<i<<" "<<index<<" "<<A.valuePtr()[j]<<endl;
  }
  p.close();
	return 0;
}