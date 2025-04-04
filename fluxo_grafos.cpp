#include <iostream>
#include <limits.h>
#include <queue>
#include <string.h>
#include <bits/stdc++.h>
using namespace std;

// Number of vertices in given graph
#define V 6
#define INF INT_MAX

int adjMatrix[V][V];
const char* vertexNames[V] = {"s", "v1", "v2", "v3", "v4", "t"};

////////////////////////////////algoritimo fordFulkerson///////////////////////////////////////

/* Returns true if there is a path from source 's' to sink
't' in residual graph. Also fills parent[] to store the
path */
bool bfs(int rGraph[V][V], int s, int t, int parent[])
{
	// Create a visited array and mark all vertices as not
	// visited
	bool visited[V];
	memset(visited, 0, sizeof(visited));

	// Create a queue, enqueue source vertex and mark source
	// vertex as visited
	queue<int> q;
	q.push(s);
	visited[s] = true;
	parent[s] = -1;

	// Standard BFS Loop
	while (!q.empty()) {
		int u = q.front();
		q.pop();

		for (int v = 0; v < V; v++) {
			if (visited[v] == false && rGraph[u][v] > 0) {
				// If we find a connection to the sink node,
				// then there is no point in BFS anymore We
				// just have to set its parent and can return
				// true
				if (v == t) {
					parent[v] = u;
					return true;
				}
				q.push(v);
				parent[v] = u;
				visited[v] = true;
			}
		}
	}

	// We didn't reach sink in BFS starting from source, so
	// return false
	return false;
}

// Returns the maximum flow from s to t in the given graph
int fordFulkerson(int graph[V][V], int s, int t)
{
	int u, v;

	// Create a residual graph and fill the residual graph
	// with given capacities in the original graph as
	// residual capacities in residual graph
	int rGraph[V]
			[V]; // Residual graph where rGraph[i][j]
				// indicates residual capacity of edge
				// from i to j (if there is an edge. If
				// rGraph[i][j] is 0, then there is not)
	for (u = 0; u < V; u++)
		for (v = 0; v < V; v++)
			rGraph[u][v] = graph[u][v];

	int parent[V]; // This array is filled by BFS and to
				// store path

	int max_flow = 0; // There is no flow initially

	// Augment the flow while there is path from source to
	// sink
	while (bfs(rGraph, s, t, parent)) {
		// Find minimum residual capacity of the edges along
		// the path filled by BFS. Or we can say find the
		// maximum flow through the path found.
		int path_flow = INT_MAX;
		for (v = t; v != s; v = parent[v]) {
			u = parent[v];
			path_flow = min(path_flow, rGraph[u][v]);
		}

		// update residual capacities of the edges and
		// reverse edges along the path
		for (v = t; v != s; v = parent[v]) {
			u = parent[v];
			rGraph[u][v] -= path_flow;
			rGraph[v][u] += path_flow;
		}

		// Add path flow to overall flow
		max_flow += path_flow;
	}

	// Return the overall flow
	return max_flow;
}

///////////////////////// algoritimo de dinic //////////////////////////////////

// A structure to represent a edge between
// two vertex
struct Edge {
	int v; // Vertex v (or "to" vertex)
		// of a directed edge u-v. "From"
		// vertex u can be obtained using
		// index in adjacent array.

	int flow; // flow of data in edge

	int C; // capacity

	int rev; // To store index of reverse
			// edge in adjacency list so that
			// we can quickly find it.
};

// Residual Graph
class Graph {
	int V; // number of vertex
	int* level; // stores level of a node
	vector<Edge>* adj;

public:
	Graph(int V)
	{
		adj = new vector<Edge>[V];
		this->V = V;
		level = new int[V];
	}

	// add edge to the graph
	void addEdge(int u, int v, int C)
	{
		// Forward edge : 0 flow and C capacity
		Edge a{ v, 0, C, (int)adj[v].size() };

		// Back edge : 0 flow and 0 capacity
		Edge b{ u, 0, 0, (int)adj[u].size() };

		adj[u].push_back(a);
		adj[v].push_back(b); // reverse edge
	}

	bool BFS(int s, int t);
	int sendFlow(int s, int flow, int t, int ptr[]);
	int DinicMaxflow(int s, int t);
};

// Finds if more flow can be sent from s to t.
// Also assigns levels to nodes.
bool Graph::BFS(int s, int t)
{
	for (int i = 0; i < V; i++)
		level[i] = -1;

	level[s] = 0; // Level of source vertex

	// Create a queue, enqueue source vertex
	// and mark source vertex as visited here
	// level[] array works as visited array also.
	list<int> q;
	q.push_back(s);

	vector<Edge>::iterator i;
	while (!q.empty()) {
		int u = q.front();
		q.pop_front();
		for (i = adj[u].begin(); i != adj[u].end(); i++) {
			Edge& e = *i;
			if (level[e.v] < 0 && e.flow < e.C) {
				// Level of current vertex is,
				// level of parent + 1
				level[e.v] = level[u] + 1;

				q.push_back(e.v);
			}
		}
	}

	// IF we can not reach to the sink we
	// return false else true
	return level[t] < 0 ? false : true;
}

// A DFS based function to send flow after BFS has
// figured out that there is a possible flow and
// constructed levels. This function called multiple
// times for a single call of BFS.
// flow : Current flow send by parent function call
// start[] : To keep track of next edge to be explored.
//		 start[i] stores count of edges explored
//		 from i.
// u : Current vertex
// t : Sink
int Graph::sendFlow(int u, int flow, int t, int start[])
{
	// Sink reached
	if (u == t)
		return flow;

	// Traverse all adjacent edges one -by - one.
	for (; start[u] < adj[u].size(); start[u]++) {
		// Pick next edge from adjacency list of u
		Edge& e = adj[u][start[u]];

		if (level[e.v] == level[u] + 1 && e.flow < e.C) {
			// find minimum flow from u to t
			int curr_flow = min(flow, e.C - e.flow);

			int temp_flow
				= sendFlow(e.v, curr_flow, t, start);

			// flow is greater than zero
			if (temp_flow > 0) {
				// add flow to current edge
				e.flow += temp_flow;

				// subtract flow from reverse edge
				// of current edge
				adj[e.v][e.rev].flow -= temp_flow;
				return temp_flow;
			}
		}
	}

	return 0;
}

// Returns maximum flow in graph
int Graph::DinicMaxflow(int s, int t)
{
	// Corner case
	if (s == t)
		return -1;

	int total = 0; // Initialize result

	// Augment the flow while there is path
	// from source to sink
	while (BFS(s, t) == true) {
		// store how many edges are visited
		// from V { 0 to V }
		int* start = new int[V + 1]{ 0 };

		// while flow is not zero in graph from S to D
		while (int flow = sendFlow(s, INT_MAX, t, start)) {

			// Add path flow to overall flow
			total += flow;
		}
	
		// Remove allocated array
		delete[] start;
	}

	// return maximum flow
	return total;
}

//////////////////////////////////////algoritimo edmonds karp///////////////////////////////////

// Função para adicionar uma aresta ao grafo
void addEdge(int from, int to, int capacity) {
    adjMatrix[from][to] = capacity;
}

// Função de Edmonds-Karp para encontrar o fluxo máximo
int edmondsKarp(int source, int sink) {
    int rGraph[V][V];  // Grafo residual
    for (int u = 0; u < V; u++) {
        for (int v = 0; v < V; v++) {
            rGraph[u][v] = adjMatrix[u][v];  // Inicializa o grafo residual
        }
    }

    int parent[V];  // Para armazenar o caminho
    int maxFlow = 0;

    // Enquanto houver um caminho do origem ao destino
    while (bfs(rGraph, source, sink, parent)) {
        int pathFlow = INF;

        // Encontrar o fluxo máximo através do caminho encontrado
        for (int v = sink; v != source; v = parent[v]) {
            int u = parent[v];
            pathFlow = std::min(pathFlow, rGraph[u][v]);
        }

        // Atualizar as capacidades e as capacidades reversas das arestas
        for (int v = sink; v != source; v = parent[v]) {
            int u = parent[v];
            rGraph[u][v] -= pathFlow;
            rGraph[v][u] += pathFlow;
        }

        maxFlow += pathFlow;

        // Código apenas para imprimir o caminho
        int path[V], pathSize = 0;
        for (int v = sink; v != -1; v = parent[v]) {
            path[pathSize++] = v;
        }

        std::cout << "Caminho: ";
        for (int i = pathSize - 1; i >= 0; i--) {
            std::cout << vertexNames[path[i]];
            if (i > 0) std::cout << " -> ";
        }
        std::cout << ", Fluxo: " << pathFlow << std::endl;
    }

    return maxFlow;
}