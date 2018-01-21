#include "Graph.h"
#include <limits>

using namespace std;
Graph::Graph(int V)
{
	this->V = V;
	adj = new list<int>[V];
}

void Graph::addEdge(int v, int w)
{
	adj[v].push_back(w);
	// add when minimum path for undirected graph
	//adj[w].push_back(v);
}

void Graph::BFS(int s)
{
	// Mark all the vertices as not visited
	bool *visited = new bool[V];
	for (int i = 0; i < V; i++)
	{
		visited[i] = false;
	}
	visited[s] = true;
	list<int> queue;
	queue.push_back(s);
	// 'i' will be used to get all adjacent vertices of a vertex
	list<int>::iterator i;
	while (!queue.empty()) 
	{
		// Dequeue a vertex from queue and print it
		s = queue.front();
		cout << s << " ";
		queue.pop_front();
		// Get all adjacent vertices of the dequeued vertex s
		// If a adjacent has not been visited, then mark it visited
		// and enqueue it
		for (i = adj[s].begin(); i != adj[s].end(); ++i)
		{
			if (!visited[*i])
			{
				visited[*i] = true;
				queue.push_back(*i);
			}
		}

	}
}

void Graph::BFS(int s, bool* visited)
{
	// Mark all the vertices as not visited
	visited[s] = true;
	list<int> queue;
	queue.push_back(s);
	// 'i' will be used to get all adjacent vertices of a vertex
	list<int>::iterator i;
	while (!queue.empty())
	{
		// Dequeue a vertex from queue and print it
		s = queue.front();
		cout << s << " ";
		queue.pop_front();
		// Get all adjacent vertices of the dequeued vertex s
		// If a adjacent has not been visited, then mark it visited
		// and enqueue it
		for (i = adj[s].begin(); i != adj[s].end(); ++i)
		{
			if (!visited[*i])
			{
				visited[*i] = true;
				queue.push_back(*i);
			}
		}

	}
}

void Graph::ConnectedComponents()
{
	bool* visited = new bool[V];
	for (int i = 0; i < V; i++)
	{
		visited[i] = false;
	}
	for (int i = 0; i < V; i++)
	{
		if (!visited[i])
		{
			BFS(i, visited);
			cout << "\n";
		}
	}

}

int Graph::Distance(int s, int v)
{
	if (s == v)
		return 0;
	int* dist = new int[V];
	bool *visited = new bool[V];
	for (int i = 0; i < V; i++)
	{
		visited[i] = false;
	}
	visited[s] = true;
	dist[s] = 0;


	list<int> queue;
	queue.push_back(s);
	// 'i' will be used to get all adjacent vertices of a vertex
	list<int>::iterator i;
	while (!queue.empty())
	{
		// Dequeue a vertex from queue and print it
		s = queue.front();
		//cout << s << " ";
		queue.pop_front();
		// Get all adjacent vertices of the dequeued vertex s
		// If a adjacent has not been visited, then mark it visited
		// and enqueue it
		for (i = adj[s].begin(); i != adj[s].end(); ++i)
		{
			if (!visited[*i])
			{
				visited[*i] = true;
				dist[*i] = dist[s]+1;
				queue.push_back(*i);
			}
		}

	}
	return dist[v];
}

vector<int> Graph::Distance(int s)
{
	vector<int> distance(V, 0);
	vector<bool> visited(V, false);

	visited[s] = true;
	distance[s] = 0;
	list<int> queue;
	queue.push_back(s);
	list<int>::iterator i;
	while (!queue.empty())
	{
		s = queue.front();
		queue.pop_front();
		for (i = adj[s].begin(); i != adj[s].end(); ++i)
		{
			if (!visited[*i])
			{
				visited[*i] = true;
				distance[*i] = distance[s] + 1;
				queue.push_back(*i);
			}
		}
	}
	return distance;
}

void Graph::DFSUtil(int s, bool* visited)
{
	visited[s] = true;
	// Mark the current node as visited and print it
	cout << s << " ";
	// Recur for all the vertices adjacent to this vertex
	list<int>::iterator i;
	for (i = adj[s].begin(); i != adj[s].end(); ++i)
	{
		if (!visited[*i])
		{
			DFSUtil(*i, visited);
		}
	}
}

void Graph::DFS(int s)
{
	bool* visited = new bool[V];
	for (int i = 0; i < V; i++)
	{
		visited[i] = false;
	}
	DFSUtil(s, visited);
}
// A recursive function used by topologicalSort
void Graph::TopologicalSortUtil(int v, bool* visited, stack<int> &Stack)
{
	// Mark the current node as visited.
	visited[v] = true;
	// Recur for all the vertices adjacent to this vertex
	list<int>::iterator i;
	for (i = adj[v].begin(); i != adj[v].end(); ++i)
	{
		if (!visited[*i])
		{
			TopologicalSortUtil(*i, visited, Stack);
		}
	}
	// Push current vertex to stack which stores result
	Stack.push(v);
}

void Graph::TopologicalSort()
{
	stack<int> Stack;
	// Mark all the vertices as not visited
	bool* visited = new bool[V];
	for (int i = 0; i < V; i++)
		visited[i] = false;
	// Call the recursive helper function to store Topological
	// Sort starting from all vertices one by one
	for (int i = 0; i < V; i++)
	{
		if (visited[i] == false)
		{
			TopologicalSortUtil(i, visited, Stack);
		}
	}
	// Print contents of stack
	while (Stack.empty() == false)
	{
		cout << Stack.top() << " ";
		Stack.pop();
	}
}
Graph Graph::GetTranspose()
{
	Graph g(V);
	for (int v = 0; v < V; v++)
	{
		// Recur for all the vertices adjacent to this vertex
		list<int>::iterator i;
		for (i = adj[v].begin(); i != adj[v].end(); ++i)
		{
			g.adj[*i].push_back(v);
		}
	}
	return g;
}
void Graph::FillOrder(int v, bool* visited, stack<int> &Stack)
{
	// Mark the current node as visited and print it
	visited[v] = true;

	// Recur for all the vertices adjacent to this vertex
	list<int>::iterator i;
	for (i = adj[v].begin(); i != adj[v].end(); ++i)
	{
		if (!visited[*i])
		{
			FillOrder(*i, visited, Stack);
		}
	}
	// All vertices reachable from v are processed by now, push v 
	Stack.push(v);
}
// The main function that finds and prints all strongly connected 
// components
void Graph::PrintSCCs()
{
	stack<int> Stack;
	// Mark all the vertices as not visited (For first DFS)
	bool *visited = new bool[V];
	for (int i = 0; i < V; i++)
	{
		visited[i] = false;
	}
	// Fill vertices in stack according to their finishing times
	for (int i = 0; i < V; i++)
	{
		if (visited[i] == false)
		{
			FillOrder(i, visited, Stack);
		}
	}

	// Create a reversed graph
	Graph gr = GetTranspose();

	// Mark all the vertices as not visited (For second DFS)
	for (int i = 0; i < V; i++)
	{
		visited[i] = false;
	}
	// Now process all vertices in order defined by Stack
	while (Stack.empty() == false)
	{
		// Pop a vertex from stack
		int v = Stack.top();
		Stack.pop();

		// print strongly connected component of the popped vertex
		if (visited[v] == false)
		{
			gr.DFSUtil(v, visited);
			cout << endl;
		}
	}
}
void Graph::shortestPath(int src)
{
	// Create a priority queue to store vertices that
	// are being preprocessed. This is weird syntax in C++.
	// Refer below link for details of this syntax
	// http://geeksquiz.com/implement-min-heap-using-stl/
	priority_queue < pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;

	// Create a vector for distances and initialize all
	// distances as infinite (INF)
	vector<int> dist(V, INT_MAX);
	map<int, int> parent;
	parent.insert(make_pair(src, NULL));
	// Insert source vertex into priority queue and initialize its distance as 0
	pq.push(make_pair(0, src));
	dist[src] = 0;

	/* Looping till priority queue becomes empty (or all
	distances are not finalized) */
	while (!pq.empty())
	{
		// The first vertex in pair is the minimum distance
		// vertex, extract it from priority queue.
		// vertex label is stored in second of pair (it
		// has to be done this way to keep the vertices
		// sorted distance (distance must be first item
		// in pair)
		int u = pq.top().second;
		pq.pop();
		// 'i' is used to get all adjacent vertices of a vertex
		list<pair<int, int>>::iterator i;
		for (i = adj[u].begin(); i != adj[u].end(); ++i)
		{
			int v = (*i).first;
			int weight = (*i).second;
			if (dist[v] > dist[u] + weight)
			{
				dist[v] = dist[u] + weight;
				if (parent.find(v) != parent.end())
				{
					parent[v] = u;
				}
				else
				{
					parent.insert(make_pair(v, u));
				}
				pq.push(make_pair(dist[v], v));
			}
		}
	}
	cout << "Vertex Distance from Source\n";
	for (int i = 0; i < V; ++i)
	{
		cout << i << " " << dist[i] << " Path: " << i <<" ";
		int j = i;
		while (parent[j] != NULL)
		{
			j = parent[j];
			cout << j << " ";
		}
		cout << src << "\n";
	}
}