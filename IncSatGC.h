#ifndef INCSATGC_INCSATGC_H
#define INCSATGC_INCSATGC_H


//main source code file, this is where the main class is located on which all operations are performed


#include <vector>
#include <iostream>
#include <chrono> //time keeping
#include <ctime> //cputime with clock()
#include <stdexcept> //exceptions
#include <limits> //get the largest numbers for int
#include <memory> //smart pointers
#include <fstream> //write to file
#include <unordered_set> //unordered set for quick member checks
#include <filesystem> //check that external binary file exists
#include <csignal> //exit signals and handing

#include "core/Solver.h" //base glucose solver
#include "utils/System.h" //cpuTime and memUsed functions
#include "Encoder.h" //encoder used to build cardinality constraints
#include <boost/process.hpp> //to run external binary to obtain a clique
#include <boost/graph/adjacency_list.hpp> //graph we use when computing connected components
#include <boost/graph/connected_components.hpp> //connected components algorithm
#include <boost/dynamic_bitset.hpp> //dynamic size bitsets for faster logic operations

#include "Graph.h"
#include "Options.h"  //options struct for different settings of the algorithm
#include "Statistics.h" //struct to store and write statistics
#include "ExtendSolvers.h" //adapts cadical to use with cardinality encodings
#include "CadicalZykovPropagator.h" //implements specific ExternalPropagator that solves problem with callbacks


//some typedefs
using Bitset = boost::dynamic_bitset<>;
using Model = Bitset; //to store truth assignment of variables
using ColorMap = std::vector<int>; //maps to each vertex the color it has been assigned, starts at color 0, -1 if uncolored
constexpr int NoColor = -1; //already defined in other header
using ConflictList = std::vector<std::tuple<int, int, int>>; //list of triplets that are the conflicts




//helper struct to store and get indices for the variables s_ij for all i<j
//stores the upper triangular matrix in a continuous vector
//rows and columns begin at 0
struct UpperTriangle{
    UpperTriangle();
    explicit UpperTriangle(int dimension, int value = -1);
    [[nodiscard]] int index(int row, int col) const ;
    void set(int row, int col, int val);
    [[nodiscard]] int get(int row, int col) const;
    [[nodiscard]] std::pair<int, int> get_ij(int index) const;

    int dimension;
    std::vector<int> sij_indices;
    std::vector<std::pair<int, int>> index_to_ij;
};

inline int UpperTriangle::index(int row, int col) const{
    assert(row != col);
    if(row > col){std::swap(row, col);}
    return row * dimension - row * (row + 1) / 2 + col - row - 1;
}

inline void UpperTriangle::set(const int row, const int col, const int val) {
    sij_indices[index(row, col)] = val;
    index_to_ij[val] = {row, col};
}

inline int UpperTriangle::get(const int row, const int col) const{
    return sij_indices[index(row, col)];
}

inline std::pair<int, int> UpperTriangle::get_ij(int index) const{
    return index_to_ij[index];
}

//main class object which is instantiated to run the algorithm
class IncSatGC {
public:
    explicit IncSatGC(const char *filename, Options options = Options());
    explicit IncSatGC(Graph::Graph in_graph, Options options_ = Options());

    ~IncSatGC() = default; //shared_ptr of solver is automatically deleted, default destructor is sufficient

    int run();

    void preprocessing();

    int compute_chromatic_number();

private:
    //the graph we are working with (might be changed in preprocessing)
    Graph::Graph graph;
    int num_vertices;
    //objects for storing interesting stats and getting the runtime options of the algorithm
    const Options options;
    Statistics stats;
    //the pointer to the solver object we are using, and the vector of literals used in assumptions
    std::shared_ptr<NSPACE::Solver> solver;
    //variables for lower/upper bound
    int lower_bound;
    int clique_lower_bound;
    int mc_lower_bound;
    int upper_bound;
    int heuristic_bound;

    //general functions to interact with used solvers in a common interface
    void new_SAT_solver();
    void reset_SAT_solver();
    void add_vars(int num_vars);
    int get_num_vars() const;
    int get_num_clauses() const;
    bool run_solver();
    void write_cnf();
    //returns the current model after a successful call to the solver
    Model get_model() const;
    //helper function to cleanly access solution value of a variable s_ij
    bool get_model_ij(const Model &model, int i, int j) const;

    //variables used in zykov tree encodings
    Graph::NeighborList complement_graph_adjacency; //to quickly iterate over non-neighbours
    UpperTriangle sij_indices; //for the pairing variables s_ij = true <=> same(i,j)
    std::vector<int> c_indices; //for the c_j variables used to count the number of colors used
    //filestream used to write full MaxSAT encoding
    std::ofstream wcnf;
    //Encoder to generate and add the clauses for at most k constraint
    openwbo::Encoder encoder;
    vec<Lit> encoder_assumptions;
    vec<Lit> cj_literals;
    vec<Lit> tmp_clause; //temporary buffer for clauses
    //idea: some c_j can trivially be set to true, remove them to reduce size of cardinality encoding
    int num_removable_cj;
    Bitset is_removable_cj;

    //store best known coloring after preprocessing and solving new SAT instance
    ColorMap current_best_coloring;
    // data computed in preprocessing and preprocessing functions
    bool solved_in_preprocessing;
    std::vector<Graph::VertexType> clique;
    Graph::Coloring heuristic_coloring;
    //use CliSAT binary to get a large clique quickly
    std::vector<Graph::VertexType> external_get_clique(bool write_graph_get_clique);
    Graph::Permutation graph_vertex_ordering; //store applied permutation
    using SubGraph = std::map<Graph::VertexType, Bitset >; //maps vertex to adjacent vertices
    int mycielsky_extension_lb(SubGraph H) const;
    //function for graph reduction, and flag whether it was successful
    bool has_removed_vertices_in_reduction = false;
    bool reduced_graph();
    //higher level preprocessing functions
    void preprocessing_clique_bound(bool write_new_graph_and_get_clique = false);
    void preprocessing_mycielsky_bound();
    void preprocessing_reductions();
    void preprocessing_initial_coloring();
    void preprocessing_clique_ordering();


    //collected notification and print functions for lb/ub updates or other related cases
    void notify_new_bound(bool res, int num_colors);
    void notify_already_used_fewer_colors(int num_colors);
    void notify_lower_bound(int num_colors);
    void notify_clique_lb(int num_colors);
    void notify_mycielsky_lb(int num_colors);
    void notify_upper_bound(int num_colors);
    void notify_heuristic_ub(int num_colors);
    void collect_bound_information(int bound, bool res);
    void print_single_k_solved_by_upper_bound(int num_colors) const;
    void print_single_k_solved_by_lower_bound(int num_colors) const;
    void print_short_cegar_stats() const;
    void print_sat_size();

    //function to verify whether a found coloring is valid
    bool is_valid_coloring(const ColorMap &colors) const;
    bool is_valid_coloring(const Graph::Coloring &coloring) const;
    //function to obtain coloring from SAT solution of direct encoding
    ColorMap obtain_coloring_from_direct_encoding(int num_colors) const;
    //functions to obtain coloring from SAT solution of zykov encoding
    int get_num_colors_from_model() const;
    ColorMap obtain_coloring_from_model() const;
    //write optimal or best found coloring to dimacs file, undoes preprocessing steps too
    void write_optimal_coloring();

    // functions that cover the single-k, top-down and bottom-up approach both for the assignment encoding and the partial order encoding
    void build_direct_encoding(int num_colors);
    bool direct_encoding_single_k();
    int direct_encoding_top_down();
    void direct_encoding_top_down_fix_variables(int full_num_colors, int colors); //helper function in top-down solving
    int direct_encoding_bottom_up();

    //functions for using the assignment encoding for the graph coloring problem
    int do_assignment_encoding();
    void build_assignment_encoding(int num_colors);
    bool assignment_encoding_single_k(); //retuns SAT or UNSAT for single k
    int assignment_encoding_top_down();
    int assignment_encoding_bottom_up();

    //functions for using the partial order encoding for the graph coloring problem
    int do_partial_order_encoding();
    void build_partial_order_encoding(int num_colors);
    bool partial_order_encoding_single_k();
    int partial_order_encoding_top_down();
    int partial_order_encoding_bottom_up();

    //functions for full encoding MaxSAT formulation
    void write_full_maxsat_encoding();
    void write_transitivity(int &count_clauses, int i, int j, int k);
    void write_all_transitivity(int &count_clauses);
    void write_cj_definition(int &count_clauses);

    // functions to initialise var indices, and building the zykov encoding + the color counting c_j + at most k
    void initialise_variable_indices(int start_index = 0);
    inline void add_transitivity(int i, int j, int k);
    void add_all_transitivity();
    void add_cj_definition();
    void add_zykov_encoding();
    void add_at_most_k(int k);
    void add_incremental_at_most_k(int k);

    // functions that cover the single-k, top-down and bottom-up approach both for the full encoding and the cegar approach
    bool zykov_encoding_run_solver();
    bool zykov_encoding_single_k();
    int zykov_encoding_top_down();
    int zykov_encoding_bottom_up();

    //functions for using the full zykov encoding. They just call the zykov_encoding_... functions
    int do_full_encoding();
    bool full_encoding_single_k();
    int full_encoding_top_down();
    int full_encoding_bottom_up();

    //functions for cegar approach of zykov encoding. They just call the zykov_encoding_... functions
    int do_cegar_approach();
    bool cegar_approach_single_k();
    int cegar_approach_top_down();
    int cegar_approach_bottom_up();

    //this function is the actual cegar loop: solve, find conflicts, solve again until UNSAT or no more conflicts are found
    bool search_add_conflicts_and_solve();
    //different methods to find existing conflicts
    ConflictList find_conflicts() const;
    ConflictList find_conflicts_naive() const;
    Graph::NeighborList get_pairing_graph(const Model &model) const;
    ConflictList find_conflicts_sparse_triangles() const;
    ConflictList find_conflicts_all_triangles() const;
    ConflictList find_conflicts_paper() const;
    void add_conflict_clauses(const ConflictList &conflicts);


    //functions for custom propagator approach using cadical and callbacks
    // for zykov encoding
    friend class CadicalZykovPropagator;
    std::unique_ptr<CadicalZykovPropagator> zykov_propagator;

    int do_zykov_propagator();
    void init_zykov_propagator();
    bool zykov_propagator_single_k();
    int zykov_propagator_top_down();
    int zykov_propagator_bottom_up();

    // function that runs cegar algorithm in a configuration as close a as possible to original paper
    int original_paper_configuration();


    //write statistics when exiting program vie user or other interrupts
    void register_write_cleanup_on_exit() const;
public:
    void write_and_cleanup();
};

static IncSatGC* INSTANCE;
static bool CalledExitBefore = false; //to avoid infinite loop if error is called in cleanup
static void write_and_cleanup_on_exit(){
    if(CalledExitBefore){
        std::exit(1);
    }
    CalledExitBefore = true;
    if(INSTANCE != nullptr){
        INSTANCE->write_and_cleanup();
        INSTANCE = nullptr;
    }
}

#endif //INCSATGC_INCSATGC_H
