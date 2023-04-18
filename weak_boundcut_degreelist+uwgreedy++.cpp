#include<Eigen/Sparse>
#include<Eigen/Dense>
#include<iostream>
#include<fstream>
#include<ctime>
#include <cstdlib>
using namespace Eigen;
using namespace std;


struct Node {                                
	int next, prev;
};


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
	SparseMatrix<int, RowMajor>A(max_size, max_size);                       
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
	int max_degree = 0;
	for (int i = 0; i < length; i++) {
		degree[i] = A.row(i).sum();   
		max_degree = max_degree > degree[i] ? max_degree : degree[i];
		if (degree[i] < bound) {
			A_cut.push_back(i);
			cut_index[i] = 1;
		}
	}
	Node* node = new Node[length];
	int* degreelist = new int[max_degree + 1]();
	for (int i = 0; i <= max_degree; i++)
		degreelist[i] = -1;
	for (int i = 0; i < length; i++) {
		node[i].next = degreelist[degree[i]];
		node[i].prev = -1;
		if (degreelist[degree[i]] != -1)
			node[degreelist[degree[i]]].prev = i;
		degreelist[degree[i]] = i;
	}
	int iters = 1;
	int* is_dealed = new int[length]();                        
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
	double time = (double)(end - start) / CLOCKS_PER_SEC;
	cout << "weak_boundcut cost time is:" << time << "s" << endl;
	//手动建立双向链表
	int start_list;
	int* next = new int[max_degree]();
	int* prev = new int[max_degree]();
	for (int i = 0; i < floor(bound); i++)
		degreelist[i] = -1;
	for (int i = floor(bound) + 1; i <= max_degree; i++) {
		if (degreelist[i] != -1) {
			start_list = i;
			prev[start_list] = -1;
			break;
		}
	}
	int last_list = start_list;
	for (int i = start_list+1; i <= max_degree; i++) {
		if (degreelist[i]!=-1) {
			prev[i] = last_list;
			next[last_list] = i;
			last_list = i;
		}
	}
	next[last_list] = -1;
	double max_density = bound;                    //bound = (double)edge_sum / (2 * node_num);
	int best_size = node_num;
	edge_sum = edge_sum/2;                       
	while (degreelist[start_list]!=-1) {
		//int i = lists[0].next;                    
		//int k = lists[i].idx;                     
		int k = degreelist[start_list];
		if (node[k].next == -1) {
			start_list = next[start_list];
			prev[start_list] = -1;
		}
		else {
			degreelist[start_list]=node[k].next;
			node[degreelist[start_list]].prev = -1;
		}
		cut_index[k] = 1;
		node_num -= 1;
		if (node_num == 0) break;
		//ans.push_back(k);
		for (int p = A.outerIndexPtr()[k]; p<A.outerIndexPtr()[k+1]; p++) { //decrement degrees of k's neighbors
			int neighbour_index = A.innerIndexPtr()[p];
			if (cut_index[neighbour_index] == 1) continue;
			edge_sum -= 1;
			//int i = itr[j];                          
			//erasenode(j);
			if (node[neighbour_index].prev == -1)
				degreelist[degree[neighbour_index]] = node[neighbour_index].next;
			else
				node[node[neighbour_index].prev].next= node[neighbour_index].next;
			if(node[neighbour_index].next != -1)
				node[node[neighbour_index].next].prev = node[neighbour_index].prev;
			//int i1 = lists[i].prev;
			int i1 = prev[degree[neighbour_index]];            
			if (degreelist[degree[neighbour_index]] == -1) {
				//eraselist(i);
				if (i1 == -1) {
					start_list = next[degree[neighbour_index]];
				}
				else
					next[i1] = next[degree[neighbour_index]];
				if (next[degree[neighbour_index]] != -1)
					prev[next[degree[neighbour_index]]] = i1;
			}
			degree[neighbour_index]--;
			//prv[j] = nxt[j] = -1;
			if (i1 == -1) {
				int i2= start_list;
				start_list = degree[neighbour_index];
				degreelist[start_list] = neighbour_index;
				next[start_list]=i2;
				prev[start_list] = -1;
				prev[i2]= start_list;
				node[neighbour_index].next = -1;
				node[neighbour_index].prev = -1;
			}
			else if (i1 != degree[neighbour_index]) {
				int i2 = next[i1];                                      
				degreelist[degree[neighbour_index]] = neighbour_index;
				next[degree[neighbour_index]] = i2;
				prev[degree[neighbour_index]] = i1;
				next[i1] = degree[neighbour_index];
				prev[i2] = degree[neighbour_index];
				node[neighbour_index].next = -1;
				node[neighbour_index].prev = -1;
			}
			else {
				node[neighbour_index].next = degreelist[i1];
				node[neighbour_index].prev = -1;
				node[degreelist[i1]].prev = neighbour_index;
				degreelist[degree[neighbour_index]] = neighbour_index;
			}
		}
		if (max_density < (double)edge_sum / node_num) {
			best_size = node_num;
			max_density = (double)edge_sum / node_num;
		}
	}
 	clock_t end1 = clock();
	time = (double)(end1 - start) / CLOCKS_PER_SEC;
	cout << "weak_boundcut+uwgreedy++ cost time is:" << time << "s" << endl;
	cout << "max_density is:"<<max_density <<"  node number is:"<<best_size<<endl;
	delete[]cut_index;
	delete[]degree;
	cout << endl;
	return 0;
}
