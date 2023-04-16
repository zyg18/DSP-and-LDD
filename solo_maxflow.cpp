#include<Eigen/Sparse>
#include<cstdint>
#include"ortools/graph/max_flow.h"
#include<Eigen/Sparse>
#include<Eigen/Dense>
#include <algorithm>
#include<iostream>
#include<fstream>
#include<ctime>
#include <cstdlib>
using namespace Eigen;
using namespace std;
using namespace operations_research;


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

double spmatrix2_density(SparseMatrix<double, RowMajor>& A, vector<int>A_index) {                    //fastest
	int length1 = A.rows();
	int length2 = A_index.size();
	int* sp_A = new int[length1]();
	for (int i = 0; i < length2; i++)
		sp_A[A_index[i]] = 1;
	double sum = 0;
	for (int i = 0; i < length2; i++) {
		int start = A.outerIndexPtr()[A_index[i]];
		int fin = A.outerIndexPtr()[A_index[i] + 1];
		for (int k = start; k < fin; k++) {
			int ind = sp_A[A.innerIndexPtr()[k]];
			if (ind != 0)
				sum += A.valuePtr()[k];
		}
	}
	double density = sum / (2 * length2);
	return density;
}



int main(int argc, char** argv)
{
	string data_in = argv[1];
	int weighted_type = atoi(argv[2]);               
	int graph_type = atoi(argv[3]);                           
	ifstream in("../../boundcut/data/" + data_in);
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
	clock_t start, end;
	start = clock();
	double downbound = 0, upbound = 0;
	for (int i = 0; i < A.rows(); i++)
		upbound = max(upbound, A.row(i).sum());


	//maxflow
	int size = A.rows();
	vector<int64_t> start_nodes, end_nodes, capacities;
	double edge_sum = A.sum() / 2;
	if (weighted_type == 0 || weighted_type == 2) {
		for (int i = 0; i < size; i++)
			for (int j = A.outerIndexPtr()[i]; j < A.outerIndexPtr()[i + 1]; j++) {
				int num = A.innerIndexPtr()[j];
				start_nodes.push_back(i);
				end_nodes.push_back(num);
			}
	}
	else if (weighted_type == 1 || weighted_type == 3) {
		for (int i = 0; i < size; i++)
			for (int j = A.outerIndexPtr()[i]; j < A.outerIndexPtr()[i + 1]; j++) {
				int num = A.innerIndexPtr()[j];
				start_nodes.push_back(i);
				end_nodes.push_back(num);
				capacities.push_back(A.valuePtr()[j]);
			}
	}
	double ss = double(size);
	double difference = 1 / (ss * (ss - 1));
	vector<int>A_maxflow_index;
	double least_density;
	cout << "upbound is: " << upbound << "  " << "downbound is: " << downbound << "  difference is:" << difference << endl;
	SimpleMaxFlow max_flow;
	while ((upbound - downbound) > difference) {
		least_density = (downbound + upbound) / 2;
		SimpleMaxFlow max_flow;
		if (weighted_type == 0 || weighted_type == 2)
			for (int i = 0; i < start_nodes.size(); ++i)
				max_flow.AddArcWithCapacity(start_nodes[i], end_nodes[i], 1);
		else if (weighted_type == 1 || weighted_type == 3) {
			for (int i = 0; i < start_nodes.size(); ++i)
				max_flow.AddArcWithCapacity(start_nodes[i], end_nodes[i], capacities[i]);
		}
		for (int i = 0; i < size; i++) {                    //The sequence number of the source point s is size, and the sequence number of the sink point is size+1
			double degree = A.row(i).sum();
			max_flow.AddArcWithCapacity(size, i, edge_sum);
			max_flow.AddArcWithCapacity(i, size + 1, edge_sum + 2 * least_density - degree);
		}
		int status = max_flow.Solve(size, size + 1);
		if (status == SimpleMaxFlow::OPTIMAL) {
			vector<int>nodes;
			max_flow.GetSourceSideMinCut(&nodes);
			auto it = find(nodes.begin(), nodes.end(), size);
			if (it != nodes.end())
				nodes.erase(it);
			sort(nodes.begin(), nodes.end());
			if (nodes.size() == 0)
				upbound = least_density;
			else {
				downbound = least_density;
				A_maxflow_index = nodes;
			}
		}
		else
			cout << "something is wrong!" << endl;
	}
	cout << "maxflow node number is:" << A_maxflow_index.size() << endl;
	double den = spmatrix2_density(A, A_maxflow_index);
	cout << "maxflow density is:" << den << endl;
	end = clock();
	double time = (double)(end - start) / CLOCKS_PER_SEC;
	cout << "maxflow cost time is:" << time << "s" << endl;
	cout << endl;
	return 0;
}