#pragma GCC optimize ("O3,unroll-loops")
#pragma GCC target ("avx2,bmi,bmi2,lzcnt,popcnt")
#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <cassert>
 
using namespace std;

const int MAX_N = 262144;
// TODO: change these to a better estimate
const int MAX_ESTIMATED_EDGES_PER_VERTEX = 1024;
// MAX_N * MAX_ESTIMATED_EDGES_PER_VERTEX
const long MAX_ESTIMATED_HEAP_SIZE = 268435456;

struct Edge {
    int vertex;
    double weight;
};

// vector<vector<Edge>> graph(MAX_N, vector<Edge>(MAX_ESTIMATED_EDGES_PER_VERTEX));
vector<vector<Edge>> graph;
double x[MAX_N], y[MAX_N], z[MAX_N], w[MAX_N];
Edge tmp_edge;

void generate_graph(int n, int d){
    mt19937_64 rng;
    uint64_t timeSeed = chrono::high_resolution_clock::now().time_since_epoch().count();
    seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32)};
    rng.seed(ss);
    uniform_real_distribution<double> unif(0, 1);

    if(d == 0){
        for(int i = 0; i < n; i++){
            for (int j = i + 1; j < n; j++){
                tmp_edge.weight = unif(rng);

                // TODO: figure out whether or not to include this edge
                tmp_edge.vertex = j;
                graph[i].push_back(tmp_edge);
                tmp_edge.vertex = i;
                graph[j].push_back(tmp_edge);
            }
        }
        return;
    }
    if(d == 2){
        for(int i = 0; i < n; i++){
            x[i] = unif(rng);
            y[i] = unif(rng);
        }
        for(int i = 0; i < n; i++){
            for(int j = i + 1; j < n; j++){
                tmp_edge.weight = sqrt((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]));

                // TODO: figure out whether or not to include this edge
                tmp_edge.vertex = j;
                graph[i].push_back(tmp_edge);
                tmp_edge.vertex = i;
                graph[j].push_back(tmp_edge);
            }
        }
        return;
    }
    if(d == 3){
        for(int i = 0; i < n; i++){
            x[i] = unif(rng);
            y[i] = unif(rng);
            z[i] = unif(rng);
        }
        for(int i = 0; i < n; i++){
            for(int j = i + 1; j < n; j++){
                tmp_edge.weight = sqrt((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]) + (z[i] - z[j]) * (z[i] - z[j]));

                // TODO: figure out whether or not to include this edge
                tmp_edge.vertex = j;
                graph[i].push_back(tmp_edge);
                tmp_edge.vertex = i;
                graph[j].push_back(tmp_edge);
            }
        }
        return;
    }
    if(d == 4){
        for(int i = 0; i < n; i++){
            x[i] = unif(rng);
            y[i] = unif(rng);
            z[i] = unif(rng);
            w[i] = unif(rng);
        }
        for(int i = 0; i < n; i++){
            for(int j = i + 1; j < n; j++){
                tmp_edge.weight = sqrt((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]) + (z[i] - z[j]) * (z[i] - z[j]) + (w[i] - w[j]) * (w[i] - w[j]));

                // TODO: figure out whether or not to include this edge
                tmp_edge.vertex = j;
                graph[i].push_back(tmp_edge);
                tmp_edge.vertex = i;
                graph[j].push_back(tmp_edge);
            }
        }
        return;
    }
}

vector<Edge> heap;

long tmp_index;
long tmp_parent_index;
long tmp_child_index;

long parent(long i){
    return (i-1)/2;
}

long child(long i){
    return 2*i+1;
}

void heap_insert(Edge edge){
    heap.push_back(edge);

    // up heapify
    tmp_index = heap.size()-1;
    if (tmp_index == 0){
        return;
    }
    tmp_parent_index = parent(tmp_index);

    while(edge.weight < heap[tmp_parent_index].weight){
        swap(heap[tmp_index], heap[tmp_parent_index]);
        tmp_index = tmp_parent_index;
        if(tmp_index == 0){
            break;
        }
        tmp_parent_index = parent(tmp_index);
    }
}

void down_heapify(){
    tmp_index = 0;
    tmp_child_index = heap[child(tmp_index)].weight < heap[child(tmp_index)+1].weight ? child(tmp_index) : child(tmp_index)+1;
    if (tmp_child_index >= heap.size()){
        return;
    }

    while(heap[tmp_index].weight > heap[tmp_child_index].weight){
        swap(heap[tmp_index], heap[tmp_child_index]);
        tmp_index = tmp_child_index;
        if (child(tmp_index) < heap.size()){
            if (child(tmp_index)+1 >= heap.size()){
                tmp_child_index = child(tmp_index);
            }
            else{
                tmp_child_index = heap[child(tmp_index)].weight < heap[child(tmp_index)+1].weight ? child(tmp_index) : child(tmp_index)+1;
            }
        }
        else{
            break;
        }
    }
}

Edge heap_extract(){
    swap(heap[0], heap[heap.size()-1]);
    Edge ret = heap.back();
    heap.pop_back();

    down_heapify();

    return ret;
}

Edge heap_insert_extract(Edge edge){
    if(heap.size() == 0 || edge.weight < heap[0].weight){
        return edge;
    }
    
    Edge ret = heap[0];
    heap[0] = edge;

    down_heapify();

    return ret;
}

// keep flipping true/false to avoid resetting array every trial
bool explored[MAX_N];
bool explored_flipped = false;

double prim_mst(int num_points){
    double mst_weight = 0;
    int current_vertex = 0;
    explored[0] = !explored_flipped;
    
    // mst contains n-1 edges
    for(int i = 0; i < num_points-1; i++){
        // cout << "currently at vertex: " << current_vertex << " adding edge number: " << i << endl;
        // for(int j = 0; j < num_points; j++) cout << explored[j] << " "; cout << endl;

        // add current edges to heap
        for(int j = 0; j < graph[current_vertex].size()-1; j++){
            if (explored_flipped ? explored[graph[current_vertex][j].vertex] : !explored[graph[current_vertex][j].vertex]){
                // cout << "inserting: " << graph[current_vertex][j].vertex << " " << graph[current_vertex][j].weight << endl;
                heap_insert(graph[current_vertex][j]);
            }
        }
        // cout << "here" << endl;
        if (explored_flipped ? explored[graph[current_vertex].back().vertex] : !explored[graph[current_vertex].back().vertex]){
            // cout << "inserting: " << graph[current_vertex].back().vertex << " " << graph[current_vertex].back().weight << endl;
            tmp_edge = heap_insert_extract(graph[current_vertex].back());
            // cout << "extracting: " << tmp_edge.vertex << " " << tmp_edge.weight << endl;
        }
        // extract until unexplored vertex
        while(explored_flipped ? !explored[tmp_edge.vertex] : explored[tmp_edge.vertex]){
            // cout << "extracting: " << tmp_edge.vertex << " " << tmp_edge.weight << endl;
            tmp_edge = heap_extract();
        }

        // update mst
        mst_weight += tmp_edge.weight;
        explored[tmp_edge.vertex] = !explored_flipped;
        current_vertex = tmp_edge.vertex;
    }

    return mst_weight;
}



int main(int argc, char* argv[]){
    // prim mst test

    // vector<Edge> v;
    // v.push_back({1, 2});
    // v.push_back({6, 4});
    // v.push_back({3, 3});
    // graph.push_back(v);
    // v.clear();

    // v.push_back({0, 2});
    // v.push_back({2, 3});
    // v.push_back({4, 2});
    // graph.push_back(v);
    // v.clear();
    
    // v.push_back({1, 3});
    // graph.push_back(v);
    // v.clear();

    // v.push_back({0, 3});
    // v.push_back({4, 5});
    // graph.push_back(v);
    // v.clear();

    // v.push_back({1, 2});
    // v.push_back({3, 5});
    // v.push_back({5, 7});
    // v.push_back({6, 6});
    // graph.push_back(v);
    // v.clear();

    // v.push_back({4, 7});
    // graph.push_back(v);
    // v.clear();

    // v.push_back({4, 6});
    // v.push_back({0, 4});
    // graph.push_back(v);

    // cout << prim_mst(7) << endl;

    int flag = atoi(argv[1]);
    int num_points = atoi(argv[2]);
    int num_trials = atoi(argv[3]);
    int dimension = atoi(argv[4]);

    double mst_weight_total = 0;
    graph.resize(num_points);

    for(int trial = 0; trial < num_trials; trial++){
        generate_graph(num_points, dimension);

        mst_weight_total += prim_mst(num_points);

        graph.clear();
        graph.resize(num_points);
        heap.clear();
        explored_flipped = !explored_flipped;
    }

    double average = mst_weight_total / num_trials;

    cout << average << " " << num_points << " " << num_trials << " " << dimension << endl;
}