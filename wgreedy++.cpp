// Densest subgraph for weighted graphs
// Uses BBST to store intermediate degrees
// 6-10 x slower than using linked lists, but simpler to understand

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include<cassert>
#include <vector>   
#include <queue> 
#include<list>
#include <set>
#include<cstring>
#include<ctime>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include<Eigen/Sparse>
#include<Eigen/Dense>
using namespace std;
using namespace Eigen;

////////////////////////////////////////////////////////////////////////////////////////
// Helper for fast input

inline char GET_CHAR() {
    const int maxn = 131072;
    static char buf[maxn], * p1 = buf, * p2 = buf;
    return p1 == p2 && (p2 = (p1 = buf) + fread(buf, 1, maxn, stdin), p1 == p2) ? EOF : *p1++;
}
inline int getInt() {
    int res(0);
    char c = GET_CHAR();
    while (c < '0') c = GET_CHAR();
    while (c >= '0') {
        res = res * 10 + (c - '0');
        c = GET_CHAR();
    }
    return res;
}

////////////////////////////////////////////////////////////////////////////////////////

vector<int> deg;

struct classcomp {
    bool operator() (const int& lhs, const int& rhs) const
    {
        return deg[lhs] < deg[rhs] || (deg[lhs] == deg[rhs] && lhs < rhs);
    }
};

set<int, classcomp> deg_sorted; //BBST storing degrees

//////////////////////////////////////////////////////////////////////////////////////////////////////
// MAIN
//////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) {

    cout << "Finding maximum subgraph density (naive version using BBST)..." << endl;

    clock_t init_begin, init_end;
    init_begin = clock();
    int iters = 1;
    string data_in = argv[1];
    int weighted_type = atoi(argv[2]);               //0表示无权图，1表示带权图，2表示处理为无权图，3表示处理为带权图
    int graph_type = atoi(argv[3]);                           //0:单部图,1：二部图
    ifstream in("./data/" + data_in);
    cout << "the data name is:" << data_in << endl;
    if (!in) {
        cout << "fail to open the file" << endl;
        return -1;
    }
    else
        cout << "opened the file" << endl;
    string comment;
    getline(in, comment);									//矩阵的解释
    char s;
    in >> s;
    int edge_num, m, n;									    //矩阵的形状和边数
    in >> edge_num >> m >> n;
    cout << "edge_num=" << edge_num << "  " << "m=" << m << "  " << "n=" << n << "  " << endl;
    int max_m = 0, max_n = 0;
    int max_size = 0;
    vector<Triplet<int>>triplets;
    if (graph_type == 0) {
        if (weighted_type == 0) {
            int a, b;
            while (in >> a >> b)                  //无权图
            {
                if (a == b)
                    continue;
                if (a < b)
                    triplets.emplace_back(a - 1, b - 1, 1);
                else
                    triplets.emplace_back(b - 1, a - 1, 1);
            }
        }
        else if (weighted_type == 1) {
            int a, b, c;
            while (in >> a >> b >> c)                  //带权图
            {
                if (a == b)
                    continue;
                if (a < b)
                    triplets.emplace_back(a - 1, b - 1, c);
                else
                    triplets.emplace_back(b - 1, a - 1, c);
            }
        }
        else if (weighted_type == 2) {          //处理为无权图
            int a, b, c;
            while (in >> a >> b >> c)
            {
                if (a == b)
                    continue;
                if (a < b)
                    triplets.emplace_back(a - 1, b - 1, 1);
                else
                    triplets.emplace_back(b - 1, a - 1, 1);
            }
        }
        else if (weighted_type == 3) {          //处理为带权图
            int a, b, c, d;
            while (in >> a >> b >> c >> d)
            {
                if (a == b)
                    continue;
                if (a < b)
                    triplets.emplace_back(a - 1, b - 1, c);
                else
                    triplets.emplace_back(b - 1, a - 1, c);
            }
        }
        else {
            cout << "Unrecognized type" << endl;
            return 0;
        }
        max_size = m > n ? m : n;
    }
    else if (graph_type == 1) {
        if (weighted_type == 0) {
            int a, b;
            while (in >> a >> b)                                           //无权图
            {
                b = m + b;
                triplets.emplace_back(a - 1, b - 1, 1);
            }
        }
        else if (weighted_type == 1) {
            int a, b, c;
            while (in >> a >> b >> c)                  //带权图
            {
                b = m + b;
                triplets.emplace_back(a - 1, b - 1, c);
            }
        }
        else if (weighted_type == 2) {          //处理为无权图
            int a, b, c;
            while (in >> a >> b >> c)
            {
                b = m + b;
                triplets.emplace_back(a - 1, b - 1, 1);
            }
        }
        else if (weighted_type == 3) {          //处理为带权图
            int a, b, c, d;
            while (in >> a >> b >> c >> d)
            {
                b = m + b;
                triplets.emplace_back(a - 1, b - 1, c);
            }
        }
        else {
            cout << "Unrecognized type" << endl;
            return 0;
        }
        max_size = m + n;
    }
    int* init_deg = new int[max_size];
    memset(init_deg, 0, sizeof(int) * max_size);
    int* w = new int[max_size];
    memset(w, 0, sizeof(int) * max_size);
    deg.resize(max_size);
    set<pair<int, int>>* nbrs = new set<pair<int, int>>[max_size];
    int sum_wts = 0;
    SparseMatrix<double, RowMajor>A(max_size, max_size);                       //构造稀疏矩阵A
    A.setFromTriplets(triplets.begin(), triplets.end());         //注意在edgelist中（4,35）有边，但是稀疏矩阵是以0为下标，所以是A(3,34)有值为1.
    A.makeCompressed();
    if (weighted_type == 0 || weighted_type == 2) {
        for (int i = 0; i < A.nonZeros(); i++)
            A.valuePtr()[i] = 1;
    }
    
    for (int i = 0; i < max_size; i++) {
        int start = A.outerIndexPtr()[i];
        int fin = A.outerIndexPtr()[i + 1];
        for (int j = start; j < fin; j++) {
            int k = A.innerIndexPtr()[j];                            //邻居的序号
            int wt = A.valuePtr()[j];
            nbrs[i].insert(make_pair(k, wt));
            nbrs[k].insert(make_pair(i, wt));
            init_deg[i] += wt;
            init_deg[k] += wt;
            sum_wts += wt;
        }
    }
   

   
    double mm_density = 0;

    vector<bool> exists(max_size);
    vector<bool> iterans(max_size);
    vector<bool> ans(max_size);

    init_end = clock();
    double init_time;
    init_time = (double)(init_end - init_begin) / CLOCKS_PER_SEC;
    cout<<sum_wts<<endl;
    cout << "Time for reading input and initialization: " << init_time << "s" << endl;

    double sum_iter_times = 0;
    clock_t run_begin, run_end;
    for (int tt = 0; tt < iters; tt++) {

        run_begin = clock();

        deg_sorted.clear();
        int deg_sum=0;
        for (int i = 0; i < max_size; i++) {
            deg[i] = w[i] + init_deg[i]; //degree for this iteration is "vertex weight" + actual degree
            deg_sorted.insert(i);
            deg_sum+=deg[i];
        }
        double iter_max_density = (double)sum_wts / max_size;
        int cur_sum_wts = sum_wts, cur_n = max_size;

        fill(exists.begin(), exists.end(), true);           //都初始化为True
        iterans = exists;                                   //都初始化为True
        int best_size = max_size;
        while (cur_n > 0) {
            cur_n--;
            int k = *(deg_sorted.begin()); //k = min degree vertex
            w[k] = deg[k]; //increment vertex weight for the next iteration (self loops)
            deg_sorted.erase(k); //delete k
            int sum_wts_decrease=cur_sum_wts;
            for (pair<int, int> j : nbrs[k]) { //decrement degrees of k's neighbors
                int nbr = j.first;
                int nbrwt = j.second;
                if (exists[nbr]) {
                    deg_sorted.erase(nbr);
                    deg[nbr] -= nbrwt;
                    cur_sum_wts -= nbrwt;
                    deg_sorted.insert(nbr);
                }
                exists[k] = false;
            }
            if (iter_max_density < (double)cur_sum_wts / cur_n) {
                iter_max_density = (double)cur_sum_wts / cur_n;
                iterans = exists;
                best_size = cur_n;
            }
        }

        if (iter_max_density > mm_density) {
            mm_density = iter_max_density;
            cout<<mm_density<<endl;
            ans = iterans;
        }

        double run_time;
        run_end = clock();
        run_time = (double)(run_end - run_begin) / CLOCKS_PER_SEC;
        sum_iter_times += run_time;
        cout << "Best size is:" << best_size << endl;
        cout << "Max density AT iteration " << tt + 1 << ": " << iter_max_density << endl;
        cout << "Max density until iteration " << tt + 1 << ": " << mm_density << endl;
        cout << "Avg time per iteration: " << sum_iter_times / (tt + 1) << " s" << endl;
        cout << "Total time: " << sum_iter_times + init_time << " s" << endl;
    }

    string output_file;
    ofstream outfile;
    return 0;
}