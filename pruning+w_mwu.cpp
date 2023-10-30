#include<Eigen/Sparse>
#include<Eigen/Dense>
#include<iostream>
#include<fstream>
#include<ctime>
#include <cstdlib>
#include<iomanip>
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
    cout<<edge_sum/2<<" "<<node_num<<endl;
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
	cout << "bound is:" << bound << " node_num is:" << node_num << endl;
	end = clock();
	double time = (double)(end - start) / CLOCKS_PER_SEC;
	cout << "weak_boundcut cost time is:" << time << "s" << endl;
	
 //node load
  A = spmatrix4_weighted(A, A_index);
	length = A.rows();
  cout<<"length is:"<<length<<endl;
	double*node_degree = new double[length];
	for(int i=0;i<length;i++)
		node_degree[i]= A.row(i).sum();
	double* node_weight = new double[length]();
  double begin_edge=A.sum()/2;
  int begin_node=length;
  //double max_density = begin_edge/begin_node;
  pair<double, int>* weight_sorted = new pair<double, int>[length];
  int T=10000;
  double edge_sum_double;
  int A_num=A.nonZeros();
  double*edge_dist=new double[A_num]();
  for(int i=0;i<A_num;i++)
    edge_dist[i]=0.5*A.valuePtr()[i];
  for(int i=0;i<length;i++)
    node_weight[i]=node_degree[i]/2;
  
  double max_density;
 	int best_num;
  ofstream p;
  p.open("./energy/mwu_"+data_in+".csv",ios::out|ios::trunc); 
  double sum_time=0;
  double E=0;
  //frank-wolfe
  for(int iter=0;iter<T;iter++){
    start=clock();
    double d;
    double*node_degree_copy=new double[length];
    double*node_weight_change=new double[length]();
    for(int i=0;i<length;i++)
      node_degree_copy[i]=node_degree[i];
  	for(int i=0;i<length; i++)
  		for (int j = A.outerIndexPtr()[i]; j < A.outerIndexPtr()[i + 1]; j++) {
  			int index = A.innerIndexPtr()[j];
        if(index<i)
          continue;
        if(node_weight[i]>node_weight[index]){
          d=1*edge_dist[j]/(iter+1);
          node_weight_change[i]-=d;
          node_weight_change[index]+=d;
          edge_dist[j]-=d;
        }
        else if(node_weight[i]<node_weight[index]){
          d=1*(A.valuePtr()[j]-edge_dist[j])/(iter+1);
          node_weight_change[i]+=d;
          node_weight_change[index]-=d;
          edge_dist[j]+=d;
        }
        else{
          d=(A.valuePtr()[j]/2-edge_dist[j])*1/(iter+1);
          node_weight_change[i]+=d;
          node_weight_change[index]-=d;
          edge_dist[j]+=d;
        }
  		}
    for(int i=0;i<length;i++)
      node_weight[i]+=node_weight_change[i];
    end=clock();
    sum_time+=(double)(end - start) / CLOCKS_PER_SEC;
    
    
    E=0;
    for(int i=0;i<length;i++)
      E+=node_weight[i]*node_weight[i];
    
    start=clock();
    //delete weighted-least node and calculate density
    max_density = begin_edge/begin_node;
    best_num=begin_node;
   	for (int i = 0; i < length; i++)
  		weight_sorted[i] = make_pair(node_weight[i], i);
      sort(weight_sorted, weight_sorted + length);
    	cut_index = new int[length]();
      edge_sum_double=begin_edge;
      node_num=begin_node;
    	for (int i = 0; i < length-1; i++){
    		int num = weight_sorted[i].second;
    		edge_sum_double -= node_degree_copy[num];
        cut_index[num]=1;
    		node_num -= 1;
    		if (edge_sum_double / node_num > max_density) {
    			max_density = edge_sum_double / node_num;
    			best_num = node_num;
    		}
    		for (int j = A.outerIndexPtr()[num]; j < A.outerIndexPtr()[num + 1]; j++) {
    			int neighbour = A.innerIndexPtr()[j];
    			if (cut_index[neighbour] == 1)
    				continue;
    			node_degree_copy[neighbour] -= A.valuePtr()[j];
    		}
    	}
   	  cout <<"iterations:"<<iter<<" node_num is:" << best_num <<" max_density is:" << max_density << endl;
      end=clock();
      time=(double)(end - start) / CLOCKS_PER_SEC;
    	cout << "sum cost time is:" << (sum_time+time) << "s" << endl; 
      p<<iter<<","<<setprecision(8)<<max_density<<","<<(sum_time+time)<<","<<best_num<<","<<E<<endl;
  }
  p.close();
	return 0;
}
