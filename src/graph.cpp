#include "graph.h"
#include <cmath>
#include <cassert>
#include <algorithm>
#include <map>
#include <unordered_set>
#include <queue>
#include <numeric> // for accumulate
#include "disjoint_set.h"
#include "heap.h"


// provide hash function for pairs
struct pair_hash{
	template <class T1, class T2>
	size_t operator () (pair<T1, T2> const &pair) const {
		size_t h1 = hash<T1>()(pair.first);
		size_t h2 = hash<T2>()(pair.second);
		return h1 ^ h2;
	}
};


// returns max{a,b}
int max(int a, int b) { return a > b ? a : b; }


// compares ((u,v),w)'s with respect to their weights
bool compare(pair<pair<int, int>, int> lhs, pair<pair<int, int>, int> rhs) {
	return lhs.second < rhs.second;
}


/////////////////////////////////////////////////////////////////////////////
//////////////////////////////CONSTRUCTORS///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


// erdos renyi constructor
Graph::Graph(int n, double p, bool i_w, bool i_d, int w_upper, int w_lower) {
	num_nodes = n;
	num_edges = 0;
	is_weighted = i_w;
	is_directed = i_d;
	vector<vector<int>> adj_list(n, vector<int>(n, -1));
	for (int i = 0; i < adj_list.size(); i++) {
		for (int j = !is_directed ? i + 1 : 0; j < adj_list.size(); j++) {
			if (i == j) {
				continue;
			}
			double r = (double)rand() / RAND_MAX; // from uniform distribution between 0 and 1
			if (r < p) { // with p probability create an undirected edge
				num_edges++;
				adj_list[i][j] = is_weighted ? rand() % w_upper + w_lower : 1;
				if (!is_directed) {
					adj_list[j][i] = adj_list[i][j];
				}
			}
		}
	}
	row_ptr = vector<int>(num_nodes+1,0);

	!is_directed ? col_ind = vector<int>(2 * num_edges) : col_ind = vector<int>(num_edges);
	weight_list = vector<int>(2*num_edges);
	int index = 0;
	for (int i = 0; i < adj_list.size(); i++) { // convert adj_list to CSR format
		for (int j = 0; j < adj_list[i].size(); j++) {
			if (adj_list[i][j] != -1) {
				row_ptr[i + 1]++;
				col_ind[index] = j;
				weight_list[index] = is_weighted ? adj_list[i][j] : 1;
				index++;
			}
		}
	}
	for (int i = 1; i < num_nodes + 1; i++) { // cumulative sum
		row_ptr[i] += row_ptr[i - 1];
	}
}


// reads edge representation file into graph
Graph::Graph(string fname) {
	ifstream input(fname.c_str());
	string line;
	// skip information lines
	do { 
		getline(input, line);
	} while (line.find("%") == !string::npos);

	istringstream ss(line);
	is_weighted = true;
	int nof_edges;
	num_edges= 0;
	ss >> num_nodes >> num_nodes >> nof_edges; // read matrix sizes, nonzero values
	vector<vector<int>> adj_list(num_nodes, vector<int>(num_nodes, -1));
	for (int i = 0; i < nof_edges; i++) {
		int v, u;
		double w;
		getline(input, line);
		ss = istringstream(line);
		ss >> v >> u;
		if (v == u) { // we do not want to have self loops
			continue;
		}
		num_edges++;
		v--;
		u--;
		adj_list[v][u] = 1;
		adj_list[u][v] = 1;
		if (is_weighted && line.find(" ") == line.rfind(" ")) { // if it has not got a weight
			is_weighted = false;
		}
		else {
			ss >> w;
			adj_list[v][u] = abs((int)w); // negative costs are not allowed
			adj_list[u][v] = abs((int)w);
		}
	}
	row_ptr = vector<int>(num_nodes + 1, 0);
	col_ind = vector<int>(2 * num_edges);
	weight_list = vector<int>(2 * num_edges);
	int index = 0;
	for (int i = 0; i < adj_list.size(); i++) { // convert adj_list into CSR format
		for (int j = 0; j < adj_list[i].size(); j++) {
			if (adj_list[i][j] != -1) {
				row_ptr[i + 1]++;
				col_ind[index] = j;
				weight_list[index] = is_weighted ? adj_list[i][j] : 1;
				index++;
			}
		}
	}
	for (int i = 1; i < num_nodes + 1; i++) { // cumulative sum
		row_ptr[i] += row_ptr[i - 1];
	}
}

/////////////////////////////////////////////////////////////////////////////
////////////////////////////UTIL FUNCTIONS///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

// adds an edge between u and v with corresponding weight
// optional: type default is undirected
void Graph::addEdge(int u, int v, int weight, string type) {
	assert(u != v);
	int m = max(u, v);
	vector<vector<int>> adj(num_nodes, vector<int>(num_nodes, -1));
	if (m + 1 > num_nodes) { // if u, or v is a new node to the graph
		adj = vector<vector<int>>(m + 1, vector<int>(m + 1, -1));
	}
	for (int k = 0; k < num_nodes; k++) { // construct an adjacency representation
		for (int edge = row_ptr[k]; edge < row_ptr[k + 1]; edge++) {
			adj[k][col_ind[edge]] = is_weighted ? weight_list[edge] : 1;
		}
	}
	adj[u][v] = weight; // mark the edge
	num_edges++;
	if (type == "undirected") { // undirected means simply two way directed edge
		adj[v][u] = weight;
		num_edges++;
	}

	if (m + 1 > num_nodes) {
		num_nodes = m + 1;
	}

	row_ptr = vector<int>(num_nodes + 1, 0);
	col_ind = vector<int>(num_edges, -1);
	weight_list = vector<int>(num_edges, -1);
	int index = 0;
	for (int i = 0; i < adj.size(); i++) { // convert adj_list to CSR format
		for (int j = 0; j < adj[i].size(); j++) {
			if (adj[i][j] != -1) {
				row_ptr[i + 1]++;
				col_ind[index] = j;
				weight_list[index] = is_weighted ? adj[i][j] : 1;
				index++;
			}
		}
	}
	for (int i = 1; i < num_nodes + 1; i++) { // cumulative sum
		row_ptr[i] += row_ptr[i - 1];
	}
}


void Graph::convertToUnWeighted() {
	is_weighted = false;
	weight_list.clear();
}


// sets graph's weights to given list
void Graph::setWeights(vector<int> & w_list) {
	assert(col_ind.size() != w_list.size());
	weight_list = w_list;
}


/////////////////////////////////////////////////////////////////////////////
/////////////////////GRAPH TRAVERSAL ALGORITHMS//////////////////////////////
/////////////////////////////////////////////////////////////////////////////


// breadth first search, returns a distance array from source node
// optional: give step size to find max step distanced nodes
vector<int> Graph::bfs(int start_node, int step_size) {
	assert(!is_weighted); // bfs works with unit weight
	vector<int> frontier(num_nodes, -1);
	vector<int> dist_arr = frontier; // every node is unvisited
	int queuestart = 0, queueend = 1, frontsize = 0;
	dist_arr[start_node] = 0; // node's distance to itself is zero
	frontier[0] = start_node; // initial frontier (queue)
	int dist = 1; // initial distance
	bool improvement = true;

	while (improvement && dist < step_size) {
		improvement = false; // if we cannot achieve to find an unvisited node we are done
		do {
			int front = frontier[queuestart++]; // choose front
			for (int edge = row_ptr[front]; edge < row_ptr[front + 1]; edge++) { // for every edge
				int adj = col_ind[edge]; // take its adjacent
				if (dist_arr[adj] == -1) { // if it is an unvisited node
					improvement = true; 
					frontier[queueend + frontsize++] = adj; // add neighbor to frontier
					dist_arr[adj] = dist; // update distance array
				}
			}
		} while (queuestart < queueend); // until current frontier finishes
		queueend += frontsize; // adjust indexes for new frontier
		frontsize = 0;
		dist++; // take a step
	}
	return dist_arr;
}


// depth first search, returns appearance order of nodes 
vector<int> Graph::dfs(int start_node) {
	vector<int> res(num_nodes, -1);
	vector<int> stack(num_nodes, -1);
	int stack_start = 0;
	stack[stack_start++] = start_node;
	int index = 1;
	res[start_node] = index;
	while (stack_start != 0) {
		int n = stack[--stack_start];
		int new_dist = res[n] + 1;
		for (int edge = row_ptr[n]; edge < row_ptr[n + 1]; edge++) {
			int adj = col_ind[edge];
			if (res[adj] == -1) {
				res[adj] = ++index;
				stack[stack_start++] = adj;
			}
		}
	}
	return res;
}

/////////////////////////////////////////////////////////////////////////////
///////////////////////////CYCLE DETECTION///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

// checks if we start from node v, we encounter a cycle
bool Graph::is_v_Cycled(int v, vector<bool> & visited, vector<bool> & stack) {
	if (visited[v] == false) {
		visited[v] = true;
		stack[v] = true;
		for (int edge = row_ptr[v]; edge < row_ptr[v + 1]; edge++) {
			int adj = col_ind[edge]; // for all neighbors of v
			if (!visited[adj] && is_v_Cycled(adj, visited, stack)) { // if neighbor is not visited and it has a cycle
				return true;
			}
			else if (stack[adj]) { // if it is already in recursive stack
				return true;
			}
		}
	}
	stack[v] = false;
	return false;
}


// checks if a graph has cycle
bool Graph::hasCycle() {
	assert(is_directed);
	vector<bool> visited(num_nodes, false); // marker for visited nodes
	vector<bool> stack(num_nodes, false); // recursion stack
	for (int v = 0; v < num_nodes; v++) {
		if (is_v_Cycled(v, visited, stack)) {
			return true;
		}
	}
	return false;
}


/////////////////////////////////////////////////////////////////////////////
//////////////////////////MINIMUM SPANNING TREE//////////////////////////////
/////////////////////////////////////////////////////////////////////////////


// returns minimum spanning tree of the graph
Graph Graph::MST() {
	DisjointSet ds; // disjoint set of number of nodes
	ds.makeSet(num_nodes);
	Graph g;
	vector<pair<pair<int, int>, int>> edges; // edge, weight
	for (int v = 0; v < num_nodes; v++) { // gather all edges
		for (int edge = row_ptr[v]; edge < row_ptr[v + 1]; edge++) {
			int adj = col_ind[edge];
			int weight = weight_list[edge];
			edges.push_back(make_pair(make_pair(v, adj), weight));
		}
	}
	sort(edges.begin(), edges.end(), compare); // ascending sort
	int new_edges = 0;
	int e = num_nodes - 1;
	for (auto x : edges) { // for all edges (u,v) in ascending order
		int a = ds.Find(x.first.first);
		int b = ds.Find(x.first.second);
		if (a != b) { // if u and v does not belong to same sets, their edge do not create a cycle
			g.addEdge(x.first.first, x.first.second, x.second);
			ds.Union(a, b);
			new_edges++;
		}
		if (new_edges == e) { // if we include enough number of edges
			return g;
		}
	}
	return g;
}


/////////////////////////////////////////////////////////////////////////////
///////////////////////SHORTEST PATH ALGORITHMS//////////////////////////////
/////////////////////////////////////////////////////////////////////////////


// dijkstra algorithm returns the shortest distance array from start_node O(logV+E)
vector<int> Graph::dijkstra(int start_node) {
	vector<int> dist_arr(num_nodes, INT_MAX);
	priority_queue<pair<int,int>> pq; // min_heap
	pq.push(make_pair(0, start_node)); // initialize minimum distance to start_node
	dist_arr[start_node] = 0;
	while (!pq.empty()) {
		int u = pq.top().second; // always choose minimum distance
		pq.pop();
		for (int edge = row_ptr[u]; edge < row_ptr[u + 1]; edge++) { // for all edges of u
			int v = col_ind[edge]; // take its neighbor
			int weight = weight_list[edge]; // and weight
			assert(weight > 0); // dijkstra does not work if weight is negative
			if (dist_arr[v] > dist_arr[u] + weight) { // if we discovered a shorter path
				dist_arr[v] = dist_arr[u] + weight; // update path
				pq.push(make_pair(dist_arr[v], v));
			}
		}
	}
	return dist_arr;
}

// belman ford algorithm returns the shortest distance array from start node O(V*E)
vector<int> Graph::bellman_ford(int start_node) {
	vector<int> dist_arr(num_nodes, INT_MAX); // at start all nodes are unreachable
	dist_arr[start_node] = 0;
	for (int v = 0; v < num_nodes-1; v++) { // for v-1 iterations
		for (int u = 0; u < num_nodes; u++) { // look and update every edge exhaustively
			for (int edge = row_ptr[u]; edge < row_ptr[u + 1]; edge++) {
				int adj = col_ind[edge];
				int weight = weight_list[edge];
				if (dist_arr[u] != INT_MAX && dist_arr[adj] > dist_arr[u] + weight) {
					dist_arr[adj] = dist_arr[u] + weight;
				}
			}
		}
	}
	return dist_arr;
}


// floyd warshall algorithm returns the shortest path matrix (from all nodes) O(V^3)
vector<vector<int>> Graph::floyd_warshall() {
	vector<vector<int>> res(num_nodes, vector<int>(num_nodes, INT_MAX/4));
	for (int i = 0; i < num_nodes; i++) { // cost of node to itself is zero
		res[i][i] = 0;
	}
	for (int v = 0; v < num_nodes; v++) { // fill adjacency matrix
		for (int edge = row_ptr[v]; edge < row_ptr[v + 1]; edge++) {
			int adj = col_ind[edge];
			int weight = weight_list[edge];
			res[v][adj] = weight;
		}
	}

	// consider all possibilities exhaustively
	for (int i = 0; i < num_nodes; i++) {
		for (int j = 0; j < num_nodes; j++) {
			for (int k = 0; k < num_nodes; k++) {
				if (res[j][i] + res[i][k] < res[j][k]) {
					res[j][k] = res[j][i] + res[i][k];
				}
			}
		}
	}
	return res;
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////GRAPH COLORING ALGORITHMS///////////////////////////////
/////////////////////////////////////////////////////////////////////////////

bool isValid(const vector<int> & color_arr, const vector<int>& row_ptr, const vector<int> & col_ind, int num_nodes) {
	for (int v = 0; v < num_nodes; v++) { // for each node v
		for (int e = row_ptr[v]; e < row_ptr[v + 1]; e++) {
			const int & adj = col_ind[e]; // for each adjacent of v
			if (color_arr[adj] == color_arr[v]) { // if color of v equals its adjacent's color
				return false;
			}
		}
	}
	return true;
}

// greedy graph coloring algorithm given order (default order[i] = i) O(V+E)
// returns the number of colors
int Graph::graphColoring(vector<int> & order) {
	if (order.size() == 0) { // if default
		order.resize(num_nodes);
		for (int i = 0; i < num_nodes; i++) { order[i] = i; }
	}
	int nof_colors = 0;
	vector<int> forbidden_arr(num_nodes, -1); // keeps colors that are forbidden to node v
	vector<int> color_arr(num_nodes, -1); // resulting color array
	for (auto v : order) { // for each node
		for (int edge = row_ptr[v]; edge < row_ptr[v + 1]; edge++) { // for each edge of node
			int adj = col_ind[edge];
			if (color_arr[adj] != -1) { // if its neighbor is already colored
				forbidden_arr[color_arr[adj]] = v; // mark its color forbidden to v
			}
		}
		for (int c = 0; c < num_nodes; c++) { // color node v with the smallest possible color (greedy)
			if (forbidden_arr[c] != v) {
				color_arr[v] = c;
				if (nof_colors < c) {
					nof_colors = c;
				}
				break;
			}
		}
	}

	return nof_colors + 1; // since we start from 0, we make +1
}


// dynamically reorders the node sequence w.r.t. 
int Graph::incidenceColoring(int start_node) {
	Heap heap(num_nodes);// each node paired with its neighbor's max color
	heap.currentSize = num_nodes;
	heap.arr[heap.locArr[start_node]].second = 1;
	heap.percolateUp(heap.locArr[start_node]);
	vector<int> forbidden_arr(num_nodes, -1); // keeps colors that are forbidden to node v
	vector<int> color_arr(num_nodes, -1); // resulting color array

	int nofcolors = 0, node_num = 0;
	while (node_num < num_nodes) {
		auto node = heap.findMax(); // node with a maximum colored neighbor
		heap.arr[heap.locArr[node.first]].second = -500;
		heap.percolateDown(heap.locArr[node.first]);
		for (int edge = row_ptr[node.first]; edge < row_ptr[node.first + 1]; edge++) { // for each edge
			int adj = col_ind[edge]; // take a neighbor
			int index = heap.locArr[adj];
			heap.arr[index].second++; // increment its neighbors val in heap
			heap.percolateUp(index);
			if (color_arr[adj] != -1) { // if it is already colored
				forbidden_arr[color_arr[adj]] = node.first; // mark its color as forbidden to node
			}
		}
		for (int c = 0; c < num_nodes; c++) { // greedily choose the smallest color that is not forbidden
			if (forbidden_arr[c] != node.first) {
				color_arr[node.first] = c;
				if (nofcolors < c) { // update if number of colors is increased
					nofcolors = c;
				}
				break;
			}
		}
		node_num++;
	}

	return nofcolors + 1;
}


// achieves a better greedy coloring performance with O(ElogV+ VlogV) complexity (heuristic)
int Graph::saturationColoring(int start_node) {

	Heap heap(num_nodes);// each node paired with its neighbor's max color
	int index = heap.locArr[start_node];
	heap.currentSize = num_nodes;
	heap.arr[index].second = 1;
	heap.percolateUp(index);
	vector<int> forbidden_arr(num_nodes, -1); // keeps colors that are forbidden to node v
	vector<int> color_arr(num_nodes, -1); // resulting color array

	int nofcolors = 0, node_num = 0, c;
	while (node_num < num_nodes) {
		auto node = heap.findMax(); // node with a maximum colored neighbor
		heap.arr[heap.locArr[node.first]].second = -500;
		heap.percolateDown(heap.locArr[node.first]);
		vector<int> neighs(row_ptr[node.first + 1] - row_ptr[node.first]); // store node's neighbors to update
		int index = 0;
		for (int edge = row_ptr[node.first]; edge < row_ptr[node.first + 1]; edge++) { // for each edge
			int adj = col_ind[edge]; // take a neighbor
			neighs[index++] = adj; // store
			if (color_arr[adj] != -1) { // if it is already colored
				forbidden_arr[color_arr[adj]] = node.first; // mark its color as forbidden to node
			}
		}
		for (c = 0; c < num_nodes; c++) { // greedily choose the smallest color that is not forbidden
			if (forbidden_arr[c] != node.first) {
				color_arr[node.first] = c;
				if (nofcolors < c) { // update if number of colors is increased
					nofcolors = c;
				}
				break;
			}
		}
		for (auto neigh : neighs) { // update each neighbor
			int index = heap.locArr[neigh];
			if (color_arr[neigh] == -1 && heap.arr[index].second < c) { // if their neighbor color increased
				heap.arr[index].second = c;
				heap.percolateUp(index); // update
			}
		}
		node_num++;
	}

	return nofcolors + 1;
}

/////////////////////////////////////////////////////////////////////////////
////////////////////////////GRAPH METRICS////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


// returns max_possible_links_of_adjacents/actual_number_of_links for each node, paired with node id
vector<pair<int, float>> Graph::clustering_coeff() {
	vector<pair<int, float>> res(num_nodes);
	for (int v = 0; v<num_nodes; v++) { // for each node
		int degree = row_ptr[v + 1] - row_ptr[v];
		int possiblelinks = degree*(degree - 1) / 2; // maximum number of links among its adjacents
		int noflinks = 0;
		unordered_set<int> set;
		for (int edge = row_ptr[v]; edge < row_ptr[v + 1]; edge++) { // insert all of its neighbors
			int adj = col_ind[edge];
			set.insert(adj);
		}
		for (int edge = row_ptr[v]; edge < row_ptr[v + 1]; edge++) {
			const int & adj = col_ind[edge]; // for each neighbor of v
			for (int link = row_ptr[adj]; link < row_ptr[adj + 1]; link++) {
				int adj_neigh = col_ind[link];// for each neighbor of neighbor
				if (set.find(adj_neigh) != set.end()) { // count number of links
					noflinks++;
				}
			}
		}
		float coeff = noflinks != 0 ? (float)possiblelinks / noflinks : possiblelinks + 1;
		res[v] = make_pair(v, coeff);
	}
	return res;
}


// returns sum_of_distance_to_all_other_nodes/num_nodes for each node
vector<pair<int, float>> Graph::closeness_centrality() {
	vector<pair<int, float>> res(num_nodes);
	vector<int> dist_arr;
	for (int v = 0; v < num_nodes; v++) {
		dist_arr = bfs(v); // take distance array for node v
		int sum_of_dist = accumulate(dist_arr.begin(), dist_arr.end(), 0); // sum of d(v,x) for all x in the graph
		float coeff = sum_of_dist > 0 ? (float)num_nodes / sum_of_dist : 0; // if coefficient is negative(meaning that graph is not connected) assign to 0
		res[v] = make_pair(v, coeff);
	}
	return res;
}


// returns degree of each node
vector<pair<int, float>> Graph::degree_order() {
	vector<pair<int, float>> res(num_nodes);
	for (int v = 0; v<num_nodes; v++) {
		res[v] = make_pair(v, row_ptr[v + 1] - row_ptr[v]);
	}
	return res;
}


// returns degree 2 of each node
vector<pair<int, float>> Graph::degree_2_order() {
	vector<pair<int, float>> res(num_nodes);
	for (int v = 0; v<num_nodes; v++) {
		vector<int> dist_arr;
		dist_arr = bfs(v, 2); // take distance array for node v
		int count = 0, val;
		for (int i = 0; i<dist_arr.size(); i++) {
			val = dist_arr[i];
			if (val == 1 || val == 2) {
				count++;
			}
		}
		res[v] = make_pair(v, count);
	}
	return res;
}

// returns degree 3 of each node
vector<pair<int, float>> Graph::degree_3_order() {
	vector<pair<int, float>> res(num_nodes);
	for (int v = 0; v<num_nodes; v++) {
		vector<int> dist_arr;
		dist_arr = bfs(v, 3); // take distance array for node v
		int count = 0, val;
		for (int i = 0; i<dist_arr.size(); i++) {
			val = dist_arr[i];
			if (val == 1 || val == 2 || val == 3) {
				count++;
			}
		}
		res[v] = make_pair(v, count);
	}
	return res;
}


vector<pair<int, float>> Graph::page_rank(int iter, float alpha) {
	vector<pair<int, float>> res(num_nodes);
	// distribute %15 evenly
	float dist_value = (1 - alpha) / num_nodes;
	for (int i = 0; i < num_nodes; i++) {
		res[i] = make_pair(i, (float)1 / num_nodes);
	}
	// initially likelyhoods are uniformly distributed
	for (int i = 0; i < iter; i++) {
		// on each iteration (required for convergence)
		vector<pair<int, float>> copy(num_nodes, pair<int, float>(0, dist_value));
		for (int j = 0; j < num_nodes; j++) {
			copy[j].first = j;
		}
		for (int v = 0; v < num_nodes; v++) { // for each node v
											  // assign total page ranks %85
			float & pr_v = copy[v].second; // update page rank of v by looking its in-degree nodes
			for (int edge = row_ptr[v]; edge < row_ptr[v + 1]; edge++) {
				// for each in degree neighbor of v (since graph is symmetric in or out degree does not matter)
				const int & adj = col_ind[edge];
				float & pr_adj = res[adj].second;
				int degree_adj = row_ptr[adj + 1] - row_ptr[adj];
				pr_v += alpha * (pr_adj / degree_adj); // page_rank_of_v <- (page_rank_of_v + page_rank_of_neighbor/out_degree_of_neighbor)
			}
		}
		res = copy;
	}
	return res;
}


/////////////////////////////////////////////////////////////////////////////
////////////////////////////OTHER PROBLEMS///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


// returns the m1*m2, O(N^3)
vector<vector<int>> matrixMultiplication(vector<vector<int>> & m1, vector<vector<int>> & m2) {
	assert(m1[0].size() == m2.size());
	vector<vector<int>> res(m1.size(), vector<int>(m2[0].size(),0)); // axb * b*c

	for (int i = 0; i < m1.size(); i++) {
		for (int j = 0; j < m2[0].size(); j++) {
			for (int k = 0; k < m1[0].size(); k++) {
				res[i][j] += m1[i][k] * m2[k][j];
			}
		}
	}
	return res;
}


// returns the number of triangles in the graph by the formula trace((adj_mat)^3)/6, O(2*V^3)
int Graph::number_of_Triangles() {
	vector<vector<int>> mat(num_nodes, vector<int>(num_nodes));
	for (int v = 0; v < num_nodes; v++) {
		for (int edge = row_ptr[v]; edge < row_ptr[v + 1]; edge++) {
			mat[v][col_ind[edge]] = 1;
		}
	}

	auto mat_cubed = matrixMultiplication(mat, matrixMultiplication(mat, mat));
	int trace = 0;
	for (int i = 0; i < mat_cubed.size(); i++) {
		trace += mat_cubed[i][i];
	}
	return trace / 6;
}



// returns maximum value and maximum valued id
void argmax(vector<int> & vec, int & id, int & val) {
	id = -1;
	val = INT_MIN;
	for (int i = 0; i < vec.size(); i++) {
		if (val < vec[i]) {
			id = i;
			val = vec[i];
		}
	}
}


// returns diameter(approximation) of a graph
int Graph::diameter() {
	int v = rand() % num_nodes; // randomly select a node
	auto dist_arr = bfs(v); // make a bfs to see the farthest node
	int id, val;
	argmax(dist_arr, id, val); // take the farthest node
	dist_arr = bfs(id); // make another bfs from it in order to go to the boundary of a graph
	argmax(dist_arr, id, val); // take the max value
	return val;
}
