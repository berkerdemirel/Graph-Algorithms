#include "disjoint_set.h"

using namespace std;

void DisjointSet::makeSet(int N) {
	for (int i = 0; i < N; i++) // create n disjoint sets
		parent[i] = i;
}

int DisjointSet::Find(int k) { // returns the root of k
	if (parent[k] == k) // if k itself a root
		return k;
	return Find(parent[k]); // recursively find
}

void DisjointSet::Union(int a, int b) {
	// find the roots of a and b
	int x = Find(a);
	int y = Find(b);
	// make one of them's parent to another
	parent[x] = y;
}

