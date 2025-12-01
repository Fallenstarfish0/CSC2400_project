# Vertex Cover Algorithm Comparison

## Project Overview

This project implements and compares three different algorithms for solving the **Vertex Cover Problem** - a classic NP-hard problem in computer science. The program provides an interactive tool to analyze algorithm performance, solution quality, and computational complexity across different graph structures.

### What is the Vertex Cover Problem?

Given an undirected graph, find the smallest set of vertices such that every edge has at least one endpoint in the set. This problem has applications in network security, resource allocation, and optimization.

## Algorithms Implemented

### 1. **Brute Force (Exact Solution)**
- **Approach**: Tests all possible vertex subsets to find the optimal solution
- **Time Complexity**: O(2^V × V × E)
- **Guarantees**: Always finds the minimum vertex cover
- **Best for**: Small graphs (≤15 vertices)

### 2. **Greedy Approximation**
- **Approach**: Repeatedly selects the vertex with the highest degree
- **Time Complexity**: O(V × E)
- **Approximation Ratio**: O(log V)
- **Best for**: Medium-sized graphs with good practical performance

### 3. **2-Approximation (Matching-Based)**
- **Approach**: Picks edges and adds both endpoints to the cover
- **Time Complexity**: O(E)
- **Approximation Ratio**: Guaranteed ≤ 2× optimal
- **Best for**: Large graphs requiring fast execution with bounded quality

## Features

- **Interactive Graph Generation**: Create random, complete, cycle, or path graphs
- **Performance Analysis**: Execution time and operation count tracking
- **Solution Quality Metrics**: Percent optimal calculations and validity checking
- **Multiple Algorithm Options**: Run all algorithms, approximations only, or fastest only
- **User-Friendly Interface**: Clear menus and comprehensive result summaries

## Graph Types Available

1. **Random Graph**: Edges added with specified probability
2. **Complete Graph**: Every vertex connected to every other vertex
3. **Cycle Graph**: Vertices connected in a ring structure
4. **Path Graph**: Vertices connected in a linear chain

## Getting Started

### Prerequisites

- **C++ Compiler**: GCC (MinGW on Windows) or any C++11 compatible compiler
- **Command Line Access**: Windows Command Prompt, PowerShell, or terminal

### Installation & Setup

1. **Clone or Download the Repository**
   ```bash
   git clone https://github.com/Fallenstarfish0/CSC2400_project
   cd CSC2400_project
   ```

2. **Navigate to Project Directory**
   ```bash
   cd "path/to/project/folder"
   ```

### Compilation

Compile the program using GCC:

```bash
g++ -o bin/project src/project.cpp
```

This creates an executable file in the `bin/` folder.

### Running the Program

Execute the compiled program:

```bash
.\bin\project
```

## Usage Guide

### 1. **Graph Type Selection**
Choose from 4 graph types:
- Enter `1` for Random Graph (specify vertices and edge probability)
- Enter `2` for Complete Graph (specify vertices)
- Enter `3` for Cycle Graph (specify vertices, minimum 3)
- Enter `4` for Path Graph (specify vertices, minimum 2)

### 2. **Algorithm Selection**
Choose which algorithms to run:
- `1`: All algorithms (recommended for ≤15 vertices)
- `2`: Approximation algorithms only (Greedy + 2-Approximation)
- `3`: 2-Approximation only (fastest option)

### 3. **Results Interpretation**

The program outputs:

#### **Individual Algorithm Results**
- Vertex cover size and actual vertices
- Execution time in milliseconds
- Basic operation count
- Optimality indication

#### **Percent Optimal Analysis**
- Solution size comparison
- Validity verification
- Percentage of optimal solution quality

#### **Performance Comparison**
- Execution time with speed ratios
- Operation count comparison
- Identification of fastest algorithm

## Example Session

```
Choose graph type:
1. Random Graph
2. Complete Graph
3. Cycle Graph
4. Path Graph
Enter choice (1-4): 1

Enter number of vertices: 10
Enter edge probability (0.0-1.0): 0.3

Choose algorithms to run:
1. All algorithms (recommended for <=15 vertices)
2. Approximation algorithms only
3. 2-Approximation only (fastest)
Enter choice (1-3): 1

============================================================
TEST: Random Graph (10 vertices, 30% edge probability)
============================================================

[Graph display and results follow...]
```

## Performance Recommendations

- **Small graphs (≤15 vertices)**: Use option 1 for complete analysis
- **Medium graphs (16-25 vertices)**: Use option 2 for approximation comparison
- **Large graphs (>25 vertices)**: Use option 3 for fastest results

## Research Applications

This tool addresses several research questions:

1. **Algorithm Performance Comparison**: How do different approximation algorithms compare in terms of solution quality and runtime?

2. **Graph Density Impact**: What is the relationship between graph density and algorithm effectiveness?

3. **Scalability Analysis**: How does performance scale with graph size across different algorithms?

## File Structure

```
CSC2400_project/
├── README.md                           # Project documentation
├── src/                                # Source code
│   └── project.cpp                     # Main C++ implementation
├── bin/                                # Compiled executables
│   └── project.exe                     # Windows executable
├── output/                             # Example program outputs
│   ├── example_out.txt                 # Sample test results
│   ├── output.txt                      # Additional test output
│   └── output_2.txt                    # Performance comparison data
├── analysis/                           # Data analysis and visualization
│   └── vertex_cover_graph_notebook.ipynb  # Jupyter analysis notebook
└── docs/                               # Academic documentation
    ├── VertexCover.pptx                # Project presentation
    └── VertexCoverReport (1).docx      # Technical report
```

## Technical Notes

- **Operation Counting**: Tracks fundamental operations (comparisons, insertions, lookups)
- **Timing Precision**: Microsecond-level timing with millisecond display
- **Memory Efficiency**: Uses optimized data structures for graph representation
- **Input Validation**: Handles invalid inputs gracefully with defaults

## Troubleshooting

### Common Issues

1. **Compilation Errors**: Ensure C++11 or later compiler
2. **Slow Performance**: Use smaller graphs or skip brute force for large inputs
3. **Character Display Issues**: Program uses ASCII characters for compatibility

### Performance Warnings

- Brute force algorithm becomes exponentially slow beyond 15 vertices
- Program warns before running potentially slow operations
- Consider using approximation-only options for large graphs

## Contributing

This project is part of CSC2400 coursework. For questions or improvements, please refer to the course guidelines.

## License

Educational use only - part of academic coursework.