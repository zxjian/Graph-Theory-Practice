#pragma once
#include <iostream>
#include <list>
#include <stack>
#include "stdafx.h"
#include <vector>
using namespace std;
// This class represents a directed graph using adjacency list representation
class Graph
{
	int V;    // No. of vertices
	list<int> *adj;    // Pointer to an array containing adjacency lists
	void DFSUtil(int v, bool* visited);  // A function used by DFS
	void TopologicalSortUtil(int v, bool* visited, stack<int> &Stack); // A function used by TopologicalSort
	// Fills Stack with vertices (in increasing order of finishing
    // times). The top element of stack has the maximum finishing 
    // time
    void FillOrder(int v, bool visited[], stack<int> &Stack);
public:
	Graph(int V);  // Constructor
	void addEdge(int v, int w); // function to add an edge to graph
	void BFS(int s);  // prints BFS traversal from a given source s
	void BFS(int s, bool* visited);
	int Distance(int s, int v);
	vector<int> Distance(int s);
	void ConnectedComponents();
	void DFS(int s); // prints DFS traversal from a given source s
	void TopologicalSort();
	// The main function that finds and prints strongly connected
	// components
	void PrintSCCs();

	// Function that returns reverse (or transpose) of this graph
	Graph GetTranspose();
	void shortestPath(int src);
};