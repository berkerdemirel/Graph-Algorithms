#ifndef DISJOINT_H
#define DISJOINT_H
#include <unordered_map>

using namespace std;

class DisjointSet {
public:
	void makeSet(int N);
	int Find(int k);
	void Union(int a, int b);
private:
	unordered_map<int, int> parent;
};
#endif