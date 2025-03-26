#ifndef IncSatGC_GRAPH_H
#define IncSatGC_GRAPH_H

/*
 * Graph.h
 * Purpose: Implementation of a class to maintain the data of a graph and calling several functions on a graph.
 * Includes reading in and building a graph as well functions to get a permuted graph, several different
 * coloring heuristics or orderings, finding a clique and functions to remove (certain) vertices from the graph.
 *
 * Modified from previous project on graph colouring
 */

#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <random>
#include <cassert>

namespace Graph{

using VertexType = int;
const VertexType UndefVertex = std::numeric_limits<VertexType>::max();
//type to store for a vertex j the set of neighboring vertices as neighbors[j]
using NeighborList = std::vector< std::vector<VertexType> >;
//note that permutations are given as positive integers, there is no i s.t. perm[i] = 0, it's always positive
using Permutation = std::vector<int>;
using ColorClass = std::set<VertexType>;
using Coloring = std::vector<ColorClass>;



//allows to print sets and vectors in a nice format to standard output
template<typename T>
std::ostream &operator<<(std::ostream &s, const std::set<T> &set);
template<typename T>
std::ostream &operator<<(std::ostream &s, const std::vector<T> &vec);

//helper function to choose a random element from a given set of elements
template<typename T>
T random_element(const std::vector<T> &v);

//returns the identity permutation, i.e. a vector with vec[i] = i of the given size
Permutation identity(int size);

//returns the inverse permutation s.t. inverse[perm[i]] = i
Permutation perm_inverse(const Permutation &perm);

//apply permutation to a set of vertices
std::vector<VertexType> apply_permutation(std::vector<VertexType> vertices, Permutation perm);


/*
 * Graph
 * Purpose: store the basic data of a graph, i.e. node count and edges and allowing to perform functions on a graph
 */
class Graph{
public:

    /*
     * Member variables
     * _ncount : the number of vertices
     * _ecount : the number of edges
     * _elist : the edges of the graph, for i = 0, ..., _ecount we have that {_elist[2i], _elist[2i+1]} is an edge
     * random_tiebreaks : if set to true, allows for randomness during choices of the coloring heuristics/orderings
     *
     * Graph : builds a graph from given vertex/edge counts and an edge list
     *         or  from a given file in either dimacs or Brendan McKay's sparse graph format ".g6" (only the first line)
     *         or  reads in a string containing a graph in Brendan McKay's sparse graph format ".g6"
     * get_neighbor_list : Puts the edge list in a more useful format, stores to a vertex
     *                     the set of adjacent vertices as neighbors[j]
     * ncount, ecount, elist : return the corresponding data field of the class
     * complement_elist : a helper function to get all the edges of the complement graph
     * get_complement : returns a graph with the same number of vertices but the edges of the complement graph
     * perm_graph : returns the permuted graph for a given permutation
     * use_random_tiebreaks : enables the use of randomness during choices in the coloring heuristics and orderings
     * print : outputs the graph in an adjacency list format
     * read_* : helper functions to build a graph from the given input source
     * simplify_graph : removes loops and duplicate edges from the graph
     *
     */

    Graph(VertexType vertex_count, VertexType edge_count, std::vector<VertexType> edge_list);
    explicit Graph(const char *filename);
    explicit Graph(std::string g6_string);

    [[nodiscard]] NeighborList get_neighbor_list() const;

    [[nodiscard]] VertexType ncount() const;
    [[nodiscard]] VertexType ecount() const;
    [[nodiscard]] double density() const;
    [[nodiscard]] std::vector<VertexType> elist() const;

    [[nodiscard]] std::vector<VertexType> complement_elist() const;
    [[nodiscard]] Graph get_complement() const;

    [[nodiscard]] Graph perm_graph(const Permutation &perm) const;

    void use_random_tiebreaks() { random_tiebreaks = true; }

    void print() const;

private:
    VertexType _ncount;
    VertexType _ecount;
    std::vector<VertexType> _elist;
    bool random_tiebreaks = false;
    //store removed vertices and vertex mapping when reducing the graph, maps v to new index remap[v]
    //careful, these are invalidated when applying graph permutations
    std::vector<VertexType> removed_vertices;
    std::vector< std::vector<VertexType> > removed_edges;
    std::vector<VertexType> remap;

    //helper functions to build the graph, either reading in dimacs or .g6 format
    void read_dimacs(const char *filename);
    void read_graph6(const char *filename);
    void read_graph6_string(std::string g6_string);

    void simplify_graph();

public:

    /*
     * Several coloring heuristics are implemented:
     *
     * dsatur: select the next vertex as the one with the highest saturation degree, ties
     *         are broken by the highest degree and randomly after that
     * dsatur_original: dsatur how it was originally introduced by Brelaz, select the next vertex
     *                  as the one with the highest saturation degree, ties are broken by the
     *                  highest degree in the uncolored subgraph and randomly after that
     * max_connected_degree: select the next vertex as the one with the highest amount of colored neighbors,
     *                       not necessarily colored differently as for dsatur.
     *                       Ties are broken by the highest degree and randomly after that.
     *                       _coloring yields a coloring, _ordering only the vertex ordering
     * constraint_graph_ordering : an ordering only algorithm, computes the ordering minimizing the
     *                             number of conflicts in the conflict graph
     *
     * Additionally, a clique can be passed to each algorithm whose vertices colors will be fixed. This can lead
     * to fewer mistakes in the heuristic since we already know these vertices will have to be colored distinctly.
     *
     * try_color_swap : implements the recolor technique for the coloring heuristics used.
     *                  if a new color class would be created, check if we can swap the colors of two vertices
     *                  to avoid having to use a new color class. This generally improves the heuristics
     *
     * Also defines OrderType, according to the algorithm by which the ordering was obtained
     */
    enum OrderType {Lexicographic, Dsatur, DsaturOriginal, MaxConnectedDegree, MinWidth};


    Coloring dsatur(Permutation &ordering, const std::vector<VertexType> &clique = {}) const;
    Coloring dsatur(Permutation &ordering, const NeighborList &neighbors, const std::vector<VertexType> &clique = {}) const;

    Coloring dsatur_original(Permutation &ordering, const std::vector<VertexType> &clique = {}) const;
    Coloring
    dsatur_original(Permutation &ordering, const NeighborList &neighbors, const std::vector<VertexType> &clique = {}) const;


    Coloring max_connected_degree_coloring(Permutation &ordering, const std::vector<VertexType> &clique = {}) const;
    Coloring max_connected_degree_coloring(Permutation &ordering, const NeighborList &neighbors,
                                           const std::vector<VertexType> &clique = {}) const;

    [[nodiscard]] Permutation max_connected_degree_ordering(const std::vector<VertexType> &clique = {}) const;
    [[nodiscard]] Permutation max_connected_degree_ordering(const NeighborList &neighbors, const std::vector<VertexType> &clique = {}) const;

    int constraint_graph_width();
    Permutation constraint_graph_ordering();

    static bool try_color_swap(VertexType max_saturated_vertex, const NeighborList &neighbors, Coloring &coloring,
                               std::vector<int> &vertex_color);

    //this is never used and not useful as a vertex ordering for this problem
    [[nodiscard]] Permutation max_degree_ordering() const;
    [[nodiscard]] Permutation max_degree_ordering(const NeighborList &neighbors) const;

    /*
     * Several functions to remove (certain) vertices from the graph.
     * note that these function are done on the graph itself and not a new one
     *
     * remove_vertex : removes a single vertex from the graph and relabels the vertices to
     *                 be labeled from 0 to n-2
     * remove_vertices : removes the given set of vertices from the graph one by one, taking care
     *                   that the right vertex is removed even after the labels change
     * remove_vertices_together : same as previous but removes the vertices all at once
     *
     * peel_graph : removes all vertices with degree strictly smaller than the specified number
     * remove_dominated_vertices : removes all dominated vertices from the graph. A vertex is dominated by another
     *                             vertex if all his neighbors are also adjacent to that other vertex
     */

    void remove_vertex(VertexType remove);
    void remove_vertices(const std::vector<VertexType> &to_remove);
    static int num_larger(VertexType v, const std::vector<VertexType>& vec);
    //just a tiny helper function to remove vertices all at once
    void remove_vertices_together(const std::vector<VertexType> &to_remove);

    int peel_graph(int peeling_degree);
    int peel_graph(int peeling_degree, const NeighborList &neighbors);

    int remove_dominated_vertices();
    int remove_dominated_vertices(const NeighborList &neighbors);


    [[nodiscard]] std::vector<VertexType> get_remapping() const;

    //recovers the reductions done on the graph by peeling and removing dominating vertices
    //Careful: only works if graph is in original ordering, i.e. no permutation was done/permutations were reverted
    std::vector<VertexType> recover_reductions();


    /*
     * function to write current graph in dimacs format to a file
     * */
    void write_dimacs(const std::string& out_path);


};

} //namespace Graph
#endif //IncSatGC_GRAPH_H
