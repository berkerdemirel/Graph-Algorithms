#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

class Graph {
public:
	Graph(bool i_w = true, bool i_d = false) :num_nodes(0), num_edges(0), is_weighted(i_w), is_directed(i_d) {}; // empty graph
	Graph(int nof_nodes, double p, bool i_w = false, bool i_d = false, int w_u=20, int w_l=1); // erdos renyi graph
	Graph(string fname); // read suite sparse graph from file
	void addEdge(int u, int v, int weight = 1, string type = "undirected");
	vector<int> bfs(int start_node, int step_size = INT_MAX); // breadth first search, returns distance array
	vector<int> dfs(int start_node); // depth first search, returns nodes appearance order
	vector<int> dijkstra(int start_node); // does not work when there is a negative cost
	vector<int> bellman_ford(int start_node); // immune to negative costs
	vector<vector<int>> floyd_warshall(); // returns pair-wise shortest distance matrix
	Graph MST(); // minimum spanning tree, using disjoint sets(kruskal algorithm)
	bool hasCycle(); // returns true if directed graph has cycle
	int number_of_Triangles(); // returns the number of triangles in the graph
	int graphColoring(vector<int> & order = vector<int>()); // colors the graph greedily given order (default, order=vertex ids)
	int incidenceColoring(int start_node = 0); // colors the graph greedily where a node is prioritized with the number of neighbors colored
	int saturationColoring(int start_node = 0); // colors the graph greedily where a node is prioritized with its neighbor's max color
	int diameter(); // returns the diameter of the graph
	vector<pair<int, float>> clustering_coeff(); // returns clustering coeff of each node
	vector<pair<int, float>> closeness_centrality(); // returns closeness centrality of each node
	vector<pair<int, float>> degree_order(); // returns degree 1 of each node
	vector<pair<int, float>> degree_2_order(); // returns degree 2 of each node
	vector<pair<int, float>> degree_3_order(); // returns degree 3 of each node
	vector<pair<int, float>> page_rank(int iter = 20, float alpha=0.85); // returns google's page rank algorithm for each node 
	bool isWeighted() { return is_weighted; }
	bool isDirected() { return is_directed; }
	void convertToUnWeighted();
	void setWeights(vector<int> & w_list);
private:
	// CSR format
	vector<int> row_ptr;
	vector<int> col_ind;
	vector<int> weight_list;
	int num_nodes;
	int num_edges;
	bool is_weighted;
	bool is_directed;
	bool is_v_Cycled(int v, vector<bool> & visited, vector<bool> & stack);
};

#endif