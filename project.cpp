#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <climits>
#include <cstdlib>
#include <ctime>

using namespace std;
using namespace chrono;

// Structure to represent an edge in the graph
struct Edge {
    int u, v;  // Two vertices that form the edge
    
    Edge(int u, int v) : u(u), v(v) {}
};

// Structure to represent a graph
class Graph {
public:
    int V;  // Number of vertices
    vector<Edge> edges;  // List of all edges
    vector<vector<int>> adjList;  // Adjacency list representation
    
    // Constructor: Initialize graph with V vertices
    Graph(int V) : V(V) {
        adjList.resize(V);
    }
    
    // Add an undirected edge between vertices u and v
    void addEdge(int u, int v) {
        edges.push_back(Edge(u, v));
        adjList[u].push_back(v);
        adjList[v].push_back(u);
    }
    
    // Print the graph structure
    void printGraph() {
        cout << "Graph with " << V << " vertices and " << edges.size() << " edges:\n";
        for (const auto& e : edges) {
            cout << "  (" << e.u << ", " << e.v << ")\n";
        }
    }
};

// Structure to store results from each algorithm
struct AlgorithmResult {
    string algorithmName;
    set<int> vertexCover;  // The vertices in the cover
    long long timeNanoseconds;  // Execution time
    bool isOptimal;  // Whether this is guaranteed to be optimal
    
    void print() const {
        cout << "\n--- " << algorithmName << " ---\n";
        cout << "Vertex Cover Size: " << vertexCover.size() << "\n";
        cout << "Vertices: { ";
        for (int v : vertexCover) {
            cout << v << " ";
        }
        cout << "}\n";
        cout << "Execution Time: " << timeNanoseconds << " nanoseconds\n";
        cout << "Optimal: " << (isOptimal ? "Yes" : "No (Approximation)") << "\n";
    }
};

// ===========================
// ALGORITHM 1: BRUTE FORCE
// ===========================
// This algorithm checks ALL possible subsets of vertices to find the minimum vertex cover
// Time Complexity: O(2^V * V * E) where V = vertices, E = edges
// Space Complexity: O(V)
// Guarantees: Always finds the OPTIMAL (minimum size) solution

class BruteForceVC {
private:
    // Helper function: Check if a given set of vertices covers all edges
    // A vertex cover is valid if every edge has at least one endpoint in the set
    bool isValidCover(const Graph& g, const set<int>& cover) {
        // Check each edge in the graph
        for (const auto& edge : g.edges) {
            // If neither endpoint of the edge is in the cover, it's not valid
            if (cover.find(edge.u) == cover.end() && 
                cover.find(edge.v) == cover.end()) {
                return false;
            }
        }
        return true;  // All edges are covered
    }
    
    // Recursive function to generate all possible subsets of vertices
    // and find the minimum valid vertex cover
    void findMinCoverRecursive(const Graph& g, int vertex, set<int>& current, 
                               set<int>& minCover, int& minSize) {
        // Pruning: If current set is already larger than best found, stop exploring
        if (current.size() >= minSize) {
            return;
        }
        
        // Base case: We've considered all vertices
        if (vertex == g.V) {
            // Check if this is a valid cover
            if (isValidCover(g, current)) {
                // Update minimum if this is smaller
                if (current.size() < minSize) {
                    minSize = current.size();
                    minCover = current;
                }
            }
            return;
        }
        
        // Recursive case 1: Include current vertex in the cover
        current.insert(vertex);
        findMinCoverRecursive(g, vertex + 1, current, minCover, minSize);
        current.erase(vertex);
        
        // Recursive case 2: Exclude current vertex from the cover
        findMinCoverRecursive(g, vertex + 1, current, minCover, minSize);
    }
    
public:
    AlgorithmResult solve(const Graph& g) {
        auto start = high_resolution_clock::now();
        
        set<int> current;  // Current subset being explored
        set<int> minCover;  // Best (minimum) cover found
        int minSize = g.V + 1;  // Initialize to impossible size
        
        // Start recursive search from vertex 0
        findMinCoverRecursive(g, 0, current, minCover, minSize);
        
        auto end = high_resolution_clock::now();
        auto duration = duration_cast<nanoseconds>(end - start).count();
        
        return {"Brute Force (Exact)", minCover, duration, true};
    }
};

// ===========================
// ALGORITHM 2: GREEDY APPROXIMATION
// ===========================
// This algorithm repeatedly selects the vertex with the highest degree (most uncovered edges)
// Time Complexity: O(V * E) where V = vertices, E = edges
// Space Complexity: O(V + E)
// Approximation Ratio: O(log V) - can be up to log(V) times larger than optimal

class GreedyVC {
public:
    AlgorithmResult solve(const Graph& g) {
        auto start = high_resolution_clock::now();
        
        set<int> cover;  // Our vertex cover solution
        set<Edge*> uncoveredEdges;  // Edges not yet covered
        
        // Initialize: All edges are uncovered
        for (const auto& edge : g.edges) {
            uncoveredEdges.insert(const_cast<Edge*>(&edge));
        }
        
        // Keep going until all edges are covered
        while (!uncoveredEdges.empty()) {
            // Count how many uncovered edges each vertex touches
            vector<int> degreeCount(g.V, 0);
            
            for (const Edge* e : uncoveredEdges) {
                degreeCount[e->u]++;
                degreeCount[e->v]++;
            }
            
            // Find vertex with maximum degree among uncovered edges
            int maxDegreeVertex = -1;
            int maxDegree = -1;
            
            for (int v = 0; v < g.V; v++) {
                if (degreeCount[v] > maxDegree) {
                    maxDegree = degreeCount[v];
                    maxDegreeVertex = v;
                }
            }
            
            // Add this vertex to our cover
            cover.insert(maxDegreeVertex);
            
            // Remove all edges covered by this vertex
            auto it = uncoveredEdges.begin();
            while (it != uncoveredEdges.end()) {
                if ((*it)->u == maxDegreeVertex || (*it)->v == maxDegreeVertex) {
                    it = uncoveredEdges.erase(it);
                } else {
                    ++it;
                }
            }
        }
        
        auto end = high_resolution_clock::now();
        auto duration = duration_cast<nanoseconds>(end - start).count();
        
        return {"Greedy Approximation", cover, duration, false};
    }
};

// ===========================
// ALGORITHM 3: 2-APPROXIMATION (MATCHING-BASED)
// ===========================
// This algorithm picks arbitrary edges and adds BOTH endpoints to the cover
// Time Complexity: O(E) where E = edges
// Space Complexity: O(V)
// Approximation Ratio: 2 - guaranteed to be at most 2x the optimal solution

class TwoApproxVC {
public:
    AlgorithmResult solve(const Graph& g) {
        auto start = high_resolution_clock::now();
        
        set<int> cover;  // Our vertex cover solution
        vector<bool> edgeCovered(g.edges.size(), false);  // Track which edges are covered
        
        // Process each edge
        for (size_t i = 0; i < g.edges.size(); i++) {
            // Skip if this edge is already covered
            if (edgeCovered[i]) continue;
            
            const Edge& e = g.edges[i];
            
            // Add both endpoints of this edge to the cover
            cover.insert(e.u);
            cover.insert(e.v);
            
            // Mark all edges incident to these vertices as covered
            for (size_t j = 0; j < g.edges.size(); j++) {
                const Edge& other = g.edges[j];
                // An edge is covered if at least one endpoint is in the cover
                if (other.u == e.u || other.u == e.v || 
                    other.v == e.u || other.v == e.v) {
                    edgeCovered[j] = true;
                }
            }
        }
        
        auto end = high_resolution_clock::now();
        auto duration = duration_cast<nanoseconds>(end - start).count();
        
        return {"2-Approximation (Matching)", cover, duration, false};
    }
};

// ===========================
// UTILITY FUNCTIONS
// ===========================

// Create a random graph with specified vertices and edge probability
Graph createRandomGraph(int V, double edgeProbability) {
    Graph g(V);
    
    // Try to add each possible edge with given probability
    for (int i = 0; i < V; i++) {
        for (int j = i + 1; j < V; j++) {
            // Random value between 0 and 1
            double randVal = (double)rand() / RAND_MAX;
            if (randVal < edgeProbability) {
                g.addEdge(i, j);
            }
        }
    }
    
    return g;
}

// Create a complete graph (every vertex connected to every other vertex)
Graph createCompleteGraph(int V) {
    Graph g(V);
    for (int i = 0; i < V; i++) {
        for (int j = i + 1; j < V; j++) {
            g.addEdge(i, j);
        }
    }
    return g;
}

// Create a cycle graph (vertices connected in a ring)
Graph createCycleGraph(int V) {
    Graph g(V);
    for (int i = 0; i < V; i++) {
        g.addEdge(i, (i + 1) % V);
    }
    return g;
}

// Create a path graph (vertices connected in a line)
Graph createPathGraph(int V) {
    Graph g(V);
    for (int i = 0; i < V - 1; i++) {
        g.addEdge(i, i + 1);
    }
    return g;
}

// Verify that a vertex cover is valid (covers all edges)
bool verifyCover(const Graph& g, const set<int>& cover) {
    for (const auto& edge : g.edges) {
        if (cover.find(edge.u) == cover.end() && 
            cover.find(edge.v) == cover.end()) {
            return false;
        }
    }
    return true;
}

// Compare and print results from all algorithms
void compareResults(const vector<AlgorithmResult>& results, const Graph& g) {
    cout << "\n" << string(60, '=') << "\n";
    cout << "COMPARISON SUMMARY\n";
    cout << string(60, '=') << "\n";
    
    // Print each algorithm's results
    for (const auto& result : results) {
        result.print();
    }
    
    // Find optimal size (from brute force)
    int optimalSize = -1;
    for (const auto& result : results) {
        if (result.isOptimal) {
            optimalSize = result.vertexCover.size();
            break;
        }
    }
    
    // Calculate approximation ratios
    cout << "\n--- Approximation Quality ---\n";
    for (const auto& result : results) {
        cout << result.algorithmName << ":\n";
        cout << "  Size: " << result.vertexCover.size();
        
        if (optimalSize > 0 && !result.isOptimal) {
            double ratio = (double)result.vertexCover.size() / optimalSize;
            cout << " (Ratio: " << fixed << setprecision(2) << ratio << "x optimal)";
        }
        cout << "\n";
        
        // Verify the cover is valid
        bool valid = verifyCover(g, result.vertexCover);
        cout << "  Valid: " << (valid ? "Yes" : "NO - ERROR!") << "\n";
    }
    
    // Compare execution times
    cout << "\n--- Performance Comparison ---\n";
    long long minTime = LLONG_MAX;
    for (const auto& result : results) {
        if (result.timeNanoseconds < minTime) {
            minTime = result.timeNanoseconds;
        }
    }
    
    for (const auto& result : results) {
        double speedup = (double)result.timeNanoseconds / minTime;
        cout << result.algorithmName << ": " 
             << result.timeNanoseconds << " ns ";
        if (speedup > 1.0) {
            cout << "(" << fixed << setprecision(2) << speedup << "x slower)";
        } else {
            cout << "(fastest)";
        }
        cout << "\n";
    }
}

// ===========================
// MAIN FUNCTION - RUN EXPERIMENTS
// ===========================

int main() {
    srand(time(0));  // Seed random number generator
    
    // Create algorithm instances
    BruteForceVC bruteForce;
    GreedyVC greedy;
    TwoApproxVC twoApprox;
    
    cout << "VERTEX COVER ALGORITHM COMPARISON\n";
    cout << "Low-Complexity Experiment: Small Graphs (≤20 vertices)\n\n";
    
    // ===========================
    // TEST 1: Small Random Graph
    // ===========================
    cout << "\n" << string(60, '=') << "\n";
    cout << "TEST 1: Random Graph (10 vertices, 30% edge probability)\n";
    cout << string(60, '=') << "\n";
    
    Graph g1 = createRandomGraph(10, 0.3);
    g1.printGraph();
    
    vector<AlgorithmResult> results1;
    results1.push_back(bruteForce.solve(g1));
    results1.push_back(greedy.solve(g1));
    results1.push_back(twoApprox.solve(g1));
    
    compareResults(results1, g1);
    
    // ===========================
    // TEST 2: Complete Graph
    // ===========================
    cout << "\n\n" << string(60, '=') << "\n";
    cout << "TEST 2: Complete Graph (8 vertices)\n";
    cout << string(60, '=') << "\n";
    
    Graph g2 = createCompleteGraph(8);
    g2.printGraph();
    
    vector<AlgorithmResult> results2;
    results2.push_back(bruteForce.solve(g2));
    results2.push_back(greedy.solve(g2));
    results2.push_back(twoApprox.solve(g2));
    
    compareResults(results2, g2);
    
    // ===========================
    // TEST 3: Cycle Graph
    // ===========================
    cout << "\n\n" << string(60, '=') << "\n";
    cout << "TEST 3: Cycle Graph (12 vertices)\n";
    cout << string(60, '=') << "\n";
    
    Graph g3 = createCycleGraph(12);
    g3.printGraph();
    
    vector<AlgorithmResult> results3;
    results3.push_back(bruteForce.solve(g3));
    results3.push_back(greedy.solve(g3));
    results3.push_back(twoApprox.solve(g3));
    
    compareResults(results3, g3);
    
    // ===========================
    // TEST 4: Path Graph
    // ===========================
    cout << "\n\n" << string(60, '=') << "\n";
    cout << "TEST 4: Path Graph (15 vertices)\n";
    cout << string(60, '=') << "\n";
    
    Graph g4 = createPathGraph(15);
    g4.printGraph();
    
    vector<AlgorithmResult> results4;
    results4.push_back(bruteForce.solve(g4));
    results4.push_back(greedy.solve(g4));
    results4.push_back(twoApprox.solve(g4));
    
    compareResults(results4, g4);
    
    // ===========================
    // TEST 5: Dense Random Graph
    // ===========================
    cout << "\n\n" << string(60, '=') << "\n";
    cout << "TEST 5: Dense Random Graph (8 vertices, 70% edge probability)\n";
    cout << string(60, '=') << "\n";
    
    Graph g5 = createRandomGraph(8, 0.7);
    g5.printGraph();
    
    vector<AlgorithmResult> results5;
    results5.push_back(bruteForce.solve(g5));
    results5.push_back(greedy.solve(g5));
    results5.push_back(twoApprox.solve(g5));
    
    compareResults(results5, g5);
    
    cout << "\n" << string(60, '=') << "\n";
    cout << "EXPERIMENT COMPLETE\n";
    cout << string(60, '=') << "\n";
    
    return 0;
}