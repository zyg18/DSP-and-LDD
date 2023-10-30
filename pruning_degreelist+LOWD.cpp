#include<Eigen/Sparse>
#include<Eigen/Dense>
#include<iostream>
#include<fstream>
#include<ctime>
#include <cstdlib>
#include<iomanip>
#include<typeinfo>
using namespace Eigen;
using namespace std;


struct Node {                                
	int next, prev;
};


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
		degree[i] = A.row(i).sum();   
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
		degreelist[deg] = -1;//Çå¿Õdegreelist[deg]
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
  int T=1000;
 	int best_num = begin_node;
  double edge_sum_double;
  int A_num=A.nonZeros();
  double*edge_dist=new double[A_num]();
  for(int i=0;i<A_num;i++)
    edge_dist[i]=0.5;
  for(int i=0;i<length;i++)
    node_weight[i]=node_degree[i]/2;
  
  double max_density;
  ofstream p;
  p.open("./output/bgreedy_"+data_in+".csv",ios::out|ios::trunc); 
  double sum_time=0;
  //BGreedy++
  for(int iter=0;iter<T;iter++){
    start=clock();
    double d;
    double*node_degree_copy=new double[length];
    for(int i=0;i<length;i++)
      node_degree_copy[i]=node_degree[i];
  	for(int i=0;i<length; i++)
  		for (int j = A.outerIndexPtr()[i]; j < A.outerIndexPtr()[i + 1]; j++) {
  			int index = A.innerIndexPtr()[j];
        if(index<i)
          continue;
        if(node_weight[i]>node_weight[index]){
          d=min((node_weight[i]-node_weight[index])/2,edge_dist[j]);
          node_weight[i]-=d;
          node_weight[index]+=d;
          edge_dist[j]-=d;
        }
        else{
          d=min((node_weight[index]-node_weight[i])/2,1-edge_dist[j]);
          node_weight[i]+=d;
          node_weight[index]-=d;
          edge_dist[j]+=d;
        }
  		}
    end=clock();
    sum_time+=(double)(end - start) / CLOCKS_PER_SEC;
    
    start=clock();
    //delete weighted-least node and calculate density
    double max_density = begin_edge/begin_node;
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
      p<<iter<<","<<setprecision(8)<<max_density<<","<<(sum_time+time)<<","<<best_num<<endl;
  }
  p.close();
	return 0;
}
