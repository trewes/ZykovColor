#ifndef INCSATGC_GRAPHMATRIX_H
#define INCSATGC_GRAPHMATRIX_H

#include <iostream>
#include <vector>
#include <set>
#include <numeric>

#include <boost/dynamic_bitset.hpp>

using Bitset = boost::dynamic_bitset<>;

#include "mnts.h"


//class that uses matrix representation for a graph
// has functionality for adding edges/contracting vertices
// and stores which vertices are in what partition
class MGraph{
public:

    MGraph() = default;
    MGraph(int n, int m, const std::vector<int> &edge_list);

    int size; //size of original graph without contractions
    int num_vertices; //keeping track of remaining vertices
    std::vector<Bitset> gmatrix; //adjacency matrix
    Bitset nodeset; //0 or 1 whether vertex is still in graph (or was contracted)
    std::vector<int> vertices; //store list of vertices for fast iteration

    //contraction data, mostly the bags for each vertex and a representation of a bag
    //we take the smallest of the vertices in a bag to be the representative
    std::vector<std::vector<int> > bag;
    std::vector<int> vertex_rep;

    //functions to check/add/remove edges. very simple, just update matrix for given u,v
    [[nodiscard]] bool has_edge(int u, int v) const;
    void add_edge(int u, int v);
    void remove_edge(int u, int v);
    //functions to check/contract vertices. adds edges to representative vertex and updates vertex bags
    [[nodiscard]] bool is_contracted(int u, int v) const;
    void contract_vertices(int u, int v);
    void separate_vertices(int u, int v);
    [[nodiscard]] double density() const;

    //data to store info for undoing edge additions and contractions done during the propagator search
    // this mainly consists of tracking the level of the graph and when which operations where performed
    int current_level = 0;
    void notify_new_level();
    // merging two vertices removes the other
    std::vector<int> removed_vertices;
    // and storing the added edges and remembering what edges where removed on what level
    std::vector<int> added_edges_trail;
    std::vector<int> added_edges_limiter;
    //and the representating vertex at each level, as well as the size of the bags to know what they like before contraction
    std::vector<std::vector<int>> vertex_rep_trail;
    std::vector<std::vector<int>> bag_sizes_trail;
    //function to restore graph data at the given level, using the added edges and partition trails
    void notify_backtrack_level(int new_level);
    //a check to make sure the graph structure makes sense
    bool check_consistency();

    //two buffers to avoid reallocation
    Bitset nv_without_nu, nu_without_nv;
    //helper function that computes the setminus of two bitsets
    static Bitset setminus(const Bitset& a, const Bitset& b);

    //prints the adjacency matrix and partition
    void print();

    //function to use sequential greedy like algorithms to compute colorings on current graph
    [[nodiscard]] int dsatur_coloring(const Bitset &clique) const;
    [[nodiscard]] int sequential_coloring(const Bitset &clique) const;
    [[nodiscard]] int IS_extract(const Bitset &clique) const;
    [[nodiscard]] int ISEQ() const;
    bool try_recolor(int vertex, std::vector<Bitset> &coloring) const;


    //function to greedily compute cliques in the graph of the current zykov node,
    // writes cliques into passed vector and returns the size of the largest clique
    // (all cliques returned are maximal and of largest size)
    int greedy_cliques(std::vector<Bitset> &clique_list, int clique_limit);
    [[nodiscard]] bool is_clique(const Bitset &clique) const;

    //struct and fuction to compute mycielsky extension of a clique in the current graph
    //helper struct used to keep track of the mycielsky subgraph being built in mycielsky_extension_clique,
    //basically a simpler version of the MGraph class itself
    struct SubGraph {
        std::vector<Bitset> matrix;
        int num_vertices = 0;
        std::vector<int> nodes;
        Bitset nodeset;
        explicit SubGraph(const int size)//creates empty subgraph with space for size nodes, i.e. empty n by n amtrix
            : matrix(size, Bitset(size)), nodeset(size) {}
        explicit SubGraph(const Bitset& clique);
    };
    //function that tries to extend a given clique into generalised mycielsky subgraph, returns number of succesful iterations
    int mycielsky_extension_clique(SubGraph &subgraph, int threshold);

    //build adjacency matrix only for the active vertices and store the mapping of indices
    std::vector<Bitset> get_active_subgraph(std::vector<int>& index_mapping) const;
    //use mnts heuristic to find clique
    int mnts_clique(std::vector<Bitset> &clique_list, int aim, int mnts_length, int mnts_depth, int random_seed = 12345) const;

};


#endif //INCSATGC_GRAPHMATRIX_H
