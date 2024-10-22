# NOI 2024

---

### Lesson 1: Voronoi Diagrams

Since Voronoi diagram construction is complex, using libraries like CGAL in C++ simplifies the process. Here’s a basic example of setting up a Voronoi diagram using the CGAL library.

```cpp
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Polygon_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K> DT;
typedef CGAL::Voronoi_diagram_2<DT> VD;

int main() {
    // Define some points
    std::vector<K::Point_2> points = {K::Point_2(0, 0), K::Point_2(1, 1), K::Point_2(2, 0)};
    DT dt(points.begin(), points.end());

    VD vd(dt);
    // Process and output Voronoi cells
    for (auto face = vd.faces_begin(); face != vd.faces_end(); ++face) {
        if (face->is_unbounded()) continue;
        std::cout << "Voronoi cell: ";
        auto v = face->dual();
        std::cout << "(" << v->point().x() << ", " << v->point().y() << ")" << std::endl;
    }
    return 0;
}
```

This code snippet demonstrates setting up Voronoi cells using CGAL. You would need to link CGAL during compilation.

---

### Lesson 2: Minimum Cut and Maximum Flow with Minimum Cost

Here's an example of implementing the Minimum Cost Maximum Flow algorithm using the Successive Shortest Path algorithm.

```cpp
#include <iostream>
#include <vector>
#include <queue>
#include <climits>

const int MAX_V = 100;
const int INF = INT_MAX;

struct Edge {
    int to, capacity, cost, rev;
};

std::vector<Edge> graph[MAX_V];
int dist[MAX_V];
int prev_v[MAX_V], prev_e[MAX_V];
int potential[MAX_V];

void add_edge(int from, int to, int capacity, int cost) {
    graph[from].push_back({to, capacity, cost, (int)graph[to].size()});
    graph[to].push_back({from, 0, -cost, (int)graph[from].size() - 1});
}

int min_cost_max_flow(int s, int t, int &flow, int V) {
    int cost = 0;
    flow = 0;
    std::fill(potential, potential + V, 0);

    while (true) {
        std::fill(dist, dist + V, INF);
        dist[s] = 0;
        std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<>> pq;
        pq.push({0, s});

        while (!pq.empty()) {
            auto [d, v] = pq.top(); pq.pop();
            if (dist[v] < d) continue;
            for (int i = 0; i < graph[v].size(); i++) {
                Edge &e = graph[v][i];
                int new_cost = dist[v] + e.cost + potential[v] - potential[e.to];
                if (e.capacity > 0 && dist[e.to] > new_cost) {
                    dist[e.to] = new_cost;
                    prev_v[e.to] = v;
                    prev_e[e.to] = i;
                    pq.push({dist[e.to], e.to});
                }
            }
        }

        if (dist[t] == INF) break;

        for (int i = 0; i < V; i++) potential[i] += dist[i];

        int aug_flow = INF;
        for (int v = t; v != s; v = prev_v[v]) {
            aug_flow = std::min(aug_flow, graph[prev_v[v]][prev_e[v]].capacity);
        }

        for (int v = t; v != s; v = prev_v[v]) {
            Edge &e = graph[prev_v[v]][prev_e[v]];
            e.capacity -= aug_flow;
            graph[v][e.rev].capacity += aug_flow;
            cost += aug_flow * e.cost;
        }
        flow += aug_flow;
    }
    return cost;
}

int main() {
    int V = 4; // number of vertices
    add_edge(0, 1, 2, 2);
    add_edge(0, 2, 1, 6);
    add_edge(1, 2, 1, 1);
    add_edge(1, 3, 1, 3);
    add_edge(2, 3, 2, 1);

    int flow;
    int min_cost = min_cost_max_flow(0, 3, flow, V);

    std::cout << "Minimum cost: " << min_cost << std::endl;
    std::cout << "Maximum flow: " << flow << std::endl;

    return 0;
}
```

This code demonstrates using the Successive Shortest Path algorithm for finding the minimum cost maximum flow in a flow network.

---

### Lesson 3: Fibonacci Heap Operations

Fibonacci heaps can be used in various algorithms, such as Dijkstra’s shortest path. Here's an example of basic Fibonacci heap operations.

```cpp
// Implementing Fibonacci Heap would require more extensive code.
// Here we provide an example of using the Boost Graph Library's Fibonacci heap.

#include <boost/heap/fibonacci_heap.hpp>
#include <iostream>

int main() {
    boost::heap::fibonacci_heap<int> fib_heap;

    fib_heap.push(10);
    fib_heap.push(20);
    fib_heap.push(30);

    std::cout << "Top element: " << fib_heap.top() << std::endl;
    fib_heap.pop();

    std::cout << "Top element after pop: " << fib_heap.top() << std::endl;
    return 0;
}
```

The Boost library provides a convenient Fibonacci heap implementation, suitable for algorithmic problems.

---

### Lesson 4: Regular Expressions and Automata

Here's a simple example for matching strings using regular expressions in C++.

```cpp
#include <iostream>
#include <regex>

int main() {
    std::string text = "hello123";
    std::regex pattern("[a-z]+[0-9]+");

    if (std::regex_match(text, pattern)) {
        std::cout << "The string matches the pattern." << std::endl;
    } else {
        std::cout << "The string does not match the pattern." << std::endl;
    }

    return 0;
}
```

This code checks if a string contains lowercase letters followed by digits.

---

### Lesson 5: Solving Systems of Linear Equations

Using Gaussian elimination to solve a system of linear equations.

```cpp
#include <iostream>
#include <vector>

using namespace std;

void gaussianElimination(vector<vector<double>>& matrix, vector<double>& result) {
    int n = matrix.size();

    for (int i = 0; i < n; i++) {
        // Partial pivoting
        for (int k = i + 1; k < n; k++) {
            if (abs(matrix[i][i]) < abs(matrix[k][i])) {
                swap(matrix[i], matrix[k]);
                swap(result[i], result[k]);
            }
        }

        // Make leading 1
        double lead = matrix[i][i];
        for (int j = 0; j < n; j++) matrix[i][j] /= lead;
        result[i] /= lead;

        // Eliminate below
        for (int k = i + 1; k < n; k++) {
            double factor = matrix[k][i];
            for (int j = 0; j < n; j++) matrix[k][j] -= factor * matrix[i][j];
            result[k] -= factor * result[i];
        }
    }

    // Back substitution
    for (int i = n - 1; i >= 0; i--) {
        for (int j = i + 1; j < n; j++) {
            result[i] -= matrix[i][j] * result[j];
        }
    }
}

int main() {
    vector<vector<double>> matrix = {
        {2, -1, 1},
        {3, 3, 9},
        {3, 3, 5}
    };
    vector<double> result = {8, -6, -4};

    gaussianElimination(matrix, result);

    for (double x : result) {
        cout << "Solution: " << x << endl;
    }
    return 0;
}
```

This code performs Gaussian elimination to solve a system of linear equations.

---
