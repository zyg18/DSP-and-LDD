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


struct Edge {                                //
    int y, next;
};

struct Node {                                //deg:点的度,next:下一个点,prev:前一个点,idx:
    int deg, next, prev, idx;
    inline void clear() {
        deg = next = prev = 0;
        idx = -1;
    }
};

Node* lists;

__inline void linklists(int x, int y) {     //点的排序(在列表里)
    if (y == 0) return;
    lists[x].next = y;
    lists[y].prev = x;
};

int* nxt, * prv, * itr;                      //用于Node类

__inline void linknodes(int x, int y) {     //点的排序
    if (y == -1) return;
    nxt[x] = y;
    prv[y] = x;
};

__inline void eraselist(int x) {            //删除点的操作（在列表里）
    lists[lists[x].prev].next = lists[x].next;
    if (lists[x].next != 0) lists[lists[x].next].prev = lists[x].prev;
};

__inline void erasenode(int x) {           //删除点的操作   erase:清除，抹去.
    if (prv[x] == -1) {
        lists[itr[x]].idx = nxt[x];
    }
    if (prv[x] != -1) nxt[prv[x]] = nxt[x];
    if (nxt[x] != -1) prv[nxt[x]] = prv[x];

};

int l = 0;
Edge* edges;
int* idx;

__inline void build(int x, int y) {         //edge[k].next记录x的下一条边,edge[k].y记录x这一次连接的对象y.
    edges[++l].next = idx[x];
    edges[l].y = y;
    idx[x] = l;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Main
int main(int argc, char** argv) {

    clock_t init_begin, init_end;
    init_begin = clock();
    int iters = 1;                 //应该是迭代次数
    string data_in = argv[1];
    int  weighted_type = atoi(argv[2]);               //0表示无权图，1表示处理为无权图
    int graph_type = atoi(argv[3]);
    ifstream in("./data/" + data_in);
    cout << "the data name is:" << data_in << endl;
    if (!in)
        cout << "fail to open the file" << endl;
    else
        cout << "opened the file" << endl;
    int edge_num, m, n;
    string comment;
    getline(in, comment);									//矩阵的解释
    char s;
    in >> s;
    in >> edge_num >> m >> n;
    edges = new Edge[edge_num * 2 + 10];
    cout << "edge_num=" << edge_num << "  " << "m=" << m << "  " << "n=" << n << "  " << endl;
    int max_m = 0, max_n = 0;
    int node_num = 0;
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
        else if (weighted_type == 1) {          //处理为无权图
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
        else {
            cout << "Unrecognized type" << endl;
            return 0;
        }
        node_num = m > n ? m : n;
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
        else if (weighted_type == 1) {          //处理为无权图
            int a, b, c;
            while (in >> a >> b >> c)
            {
                b = m + b;
                triplets.emplace_back(a - 1, b - 1, 1);
            }
        }
        else {
            cout << "Unrecognized type" << endl;
            return 0;
        }
        node_num = m + n;
    }
    idx = new int[node_num];
    memset(idx, 0, sizeof(int) * node_num);
    int* init_deg = new int[node_num];
    memset(init_deg, 0, sizeof(int) * node_num);
    l = 0;
    lists = new Node[node_num + 2 * edge_num + 10];
    int n_list = 0;
    itr = new int[node_num];
    int* deg = new int[node_num], * w = new int[node_num], * pos = new int[node_num];
    memset(deg, 0, sizeof(int) * node_num);
    memset(w, 0, sizeof(int) * node_num); //initial vertex weights=0, i.e., no self loops at the start 
    memset(pos, 0, sizeof(int) * node_num);
    prv = new int[node_num]; nxt = new int[node_num];
    SparseMatrix<double, RowMajor>A(node_num, node_num);                       //构造稀疏矩阵A
    A.setFromTriplets(triplets.begin(), triplets.end());         //注意在edgelist中（4,35）有边，但是稀疏矩阵是以0为下标，所以是A(3,34)有值为1.
    A.makeCompressed();
    for (int i = 0; i < A.nonZeros(); i++)
        A.valuePtr()[i] = 1;
    for (int i = 0; i < node_num; i++) {
        int start = A.outerIndexPtr()[i];
        int fin = A.outerIndexPtr()[i + 1];
        for (int j = start; j < fin; j++) {
            int k = A.innerIndexPtr()[j];
            build(i, k);
            build(k, i);
            init_deg[i]++;
            init_deg[k]++;
        }
    }




    pair<int, int>* deg_sorted = new pair<int, int>[node_num];
    vector<int> m_ans;
    double mm_density = 0;

    init_end = clock();
    double sum_iter_times = 0;
    double init_time;
    init_time = (double)(init_end - init_begin) / CLOCKS_PER_SEC;
    cout << "Time for reading input and initialization: " << init_time << " s" << endl;
    clock_t run_begin, run_end;
    for (int tt = 0; tt < iters; tt++) {
        run_begin = clock();
        for (int i = 0; i < node_num; i++) {
            nxt[i] = prv[i] = -1;
            pos[i] = 0;

            deg[i] = w[i] + init_deg[i];                                       //degree for this iteration is "vertex weight" + actual degree
            deg_sorted[i] = make_pair(deg[i], i);
        }
        sort(deg_sorted, deg_sorted + node_num);                                    //对点的度进行排序,对pair里的first进行升序排列。
        n_list = 0;

        for (int i = 0; i < node_num; i++) {
            int v = deg_sorted[i].second;                                     //v是度最小的点
            if (n_list == 0 || lists[n_list].deg != deg_sorted[i].first) {     //应该是让lists里依次填入最小的度
                ++n_list;
                lists[n_list].clear();
                linklists(n_list - 1, n_list);
                lists[n_list].deg = deg_sorted[i].first;
            }

            linknodes(v, lists[n_list].idx);
            lists[n_list].idx = v;
            itr[v] = n_list;
        }
        int edge_sum = A.nonZeros();
        double max_density = (double)edge_sum / node_num;
        int cur_m = edge_sum, cur_n = node_num;
        vector<int> ans;
        int max_size = 0;
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
            ans.push_back(k);
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
        reverse(ans.begin(), ans.end());
        ans.resize(node_num - max_size);
        if (max_density > mm_density) {
            m_ans = ans;
            mm_density = max_density;
        }
        double run_time;
        run_end = clock();
        run_time = (double)(run_end - run_begin) / CLOCKS_PER_SEC;
        sum_iter_times += run_time;
        cout << "Max size is:" << max_size << endl;
        cout << "Max density until iteration " << tt + 1 << ": " << mm_density << endl;
        cout << "Avg time per iteration: " << sum_iter_times / (tt + 1) << "s" << endl;
        cout << "Total time: " << sum_iter_times + init_time << " s" << endl;
    }
    return 0;
}