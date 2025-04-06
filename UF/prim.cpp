// C++ Program to implement the
// Prim's Algorithm using Adjecency matrix
#include <bits/stdc++.h>
using namespace std;

// Function to construct and print the MST
void primMST(vector<vector<int>> graph) {
    int v = graph.size();
  
    // vector to store the parent of vertex
    vector<int> parent(v);

    // vector holds the weight/ cost of the MST
    vector<int> key(v);

    // vector to represent the set of
    // vertices included in MST
    vector<bool> vis(v);

    priority_queue<pair<int, int>, 
    vector<pair<int, int>>, 
    greater<pair<int, int>>> pq;
  
    // Initialize all key vector as INFINITE
    // and vis vector as false
    for (int i = 0; i < v; i++) {
        key[i] = INT_MAX;
        vis[i] = false;
    }

    // Always include the first vertex in MST.
    // Make key 0 so that this vertex is
    // picked as the first vertex.
    key[0] = 0;

    // First node is always the root of MST
    parent[0] = -1;

    // Push the source vertex to the min-heap
    pq.push({0, 0});

    while (!pq.empty()) {
        int node = pq.top().second;
        pq.pop();
        vis[node] = true;
        for (int i = 0; i < v; i++) {
          
            // If the vertex is not visited
            // and the edge weight of neighbouring
            // vertex is less than key value of
            // neighbouring vertex then update it.
            if (!vis[i] && graph[node][i] != 0
                && graph[node][i] < key[i]) {
                pq.push({graph[node][i], i});
                key[i] = graph[node][i];
                parent[i] = node;
            }
        }
    }

    // Print the edges and their
    // weights in the MST
    cout << "Edge \tWeight\n";
    for (int i = 1; i < v; i++) {
        cout << parent[i] << " - " << i 
          << " \t" << graph[i][parent[i]] << " \n";
    }
}