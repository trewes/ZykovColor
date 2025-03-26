#ifndef INCSATGC_OPTIONS_H
#define INCSATGC_OPTIONS_H

#include <iostream>
#include <filesystem>
#include <optional> //optional values for variables
#include <boost/program_options.hpp> //use boost to parse options
namespace po = boost::program_options;

struct Options {
    //definition of all enums and corresponding names as static members of the struct
    enum Encoding {
        AssignmentEncoding, FullMaxSAT, FullEncoding, CEGAR, PartialOrderEncoding, ZykovPropagator
    };
    inline static const std::vector<std::string> tostr_encoding{
        "Assignment Encoding", "Full MaxSAT", "Full Encoding", "CEGAR", "Partial Order Encoding",
        "Zykov Encoding Propagator"
    };
    enum CheckAlgorithm {
        NaiveChecker, SparseTrianglesChecker, AllTrianglesChecker, PaperChecker
    };
    inline static const std::vector<std::string> tostr_checker{
        "Naive Checker", "Sparse Triangles Checker", "All Triangles Checker", "Paper Checker"
    };
    enum SearchStrategy {
        TopDown, BottomUp, SingleK
    };
    inline static const std::vector<std::string> tostr_strategy{
        "Top-Down", "Bottom-Up", "Single K"
    };
    enum Solver {
        Glucose, CaDiCaL
    };
    inline static const std::vector<std::string> tostr_solver{
        "Glucose", "CaDiCaL"
    };
    enum Verbosity {
        Quiet, Normal, Verbose, Debug
    };
    inline static const std::vector<std::string> tostr_verbosity{
        "Quiet", "Normal", "Verbose", "Debug"
    };

    Options(); //sets default values
    Options(int argc, char** argv); //parses command line for arguments
    void print_header() const; //print small program header with machine name and pid
    void print() const; //print what options are set

    //some information about what command is run and with what instance file
    std::string full_cmd;
    std::string filepath;
    std::string filename;
    //algorithm configuration
    Encoding encoding;
    CheckAlgorithm checker;
    SearchStrategy strategy;
    Solver solver;
    Verbosity verbosity;
    std::optional<int> specific_num_colors;
    bool disable_preprocessing;
    bool reduce_graph;
    bool use_clique_in_ordering;
    bool use_mycielsky_lb;
    bool remove_trivial_cj;
    bool assignment_encoding_amo;
    bool write_cnf_only;
    //options for zykov propagator
    enum ZykovPropagatorDecisionStrategy {
        CadicalZykov, FirstLiteral, ISUN, ImitateDsatur, BagSize
    };
    inline static const std::vector<std::string> tostr_zykov_strategy{
        "Cadical", "FirstLiteral", "ISUN", "ImitateDsatur", "BagSize"
    };
    ZykovPropagatorDecisionStrategy zykov_propagator_decision_strategy;
    bool disable_cardinality_constraints;
    enum ZykovPropagatorColoringAlgorithm {
        None, FastDsatur, SortedSEQ, IteratedIS, IteratedSEQ
    };
    inline static const std::vector<std::string> tostr_color_algorithm{
        "None", "FastDsatur", "SortedSEQ", "IteratedIS", "IteratedSEQ"
    };
    ZykovPropagatorColoringAlgorithm zykov_coloring_algorithm;
    //general propagator options, now that Assignment propagator also uses pruning as in Zykov approach
    int prop_clique_limit;
    bool use_clique_explanation_clauses;
    bool use_mycielsky_explanation_clauses;
    int mycielsky_threshold;
    bool use_dominated_vertex_decisions;
    bool enable_positive_pruning;
    bool enable_negative_pruning;
    //disable incremental bottom-up solving for Zykov encoding to test difference
    bool non_incremental_solving;
    //option that sets configuration to be as close as possible to description in original paper
    bool original_paper_configuration;
    //if given, write or append collected statistics to a csv file
    std::string stats_csvfile;
    //optional filepath of where to write a found coloring
    std::string coloringfilepath;
    //optionally enable detailed backtracking stats for propagators
    bool enable_detailed_backtracking_stats;

    //MNTS parameters
    int mnts_length;
    int mnts_depth;

private:
    static std::string enum_names_to_string(const std::vector<std::string>& enum_strings);
    static std::string option_description(std::string desc, const std::vector<std::string> &enum_strings = {});
};






#endif //INCSATGC_OPTIONS_H
