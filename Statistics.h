#ifndef INCSATGC_STATISTICS_H
#define INCSATGC_STATISTICS_H


#include <iostream>
#include <fstream> //write to file
#include <chrono> //time keeping
#include <ctime> //cputime with clock()
#include <sys/resource.h> //cputime via resource usage linux

#include "Options.h"  //options struct for different settings of the algorithm

using Duration = std::chrono::duration<double>;
using TimePoint = double; //std::clock_t


//statistics struct
class Statistics {
public:

    explicit Statistics(const Options &options);

    //variables and functions for convenient timekeeping
    enum Phase : std::uint8_t {
        Preprocessing,
            PreprocessingClique,
            PreprocessingMycielsky,
            PreprocessingReductions,
            PreprocessingInitialColoring,
        Algorithm,
            BuildEncoding,
            BuildAtMostK,
            CEGAR,
            SatSolver,
                PropagatorAssignment,
                PropagatorBacktrack,
                PropagatorDecide,
                PropagatorComputeCliques,
                PropagatorCliqueClauses,
                PropagatorMycielskyClauses,
                PropagatorPositivePruning,
                PropagatorNegativePruning,
                TestColorTime,
        Total
    };

    void start_phase(Phase phase);
    void end_phase(Statistics::Phase phase);
    [[nodiscard]] Duration duration_of(Phase phase) const;
    [[nodiscard]] double duration_in_sec(Phase phase) const;
    [[nodiscard]] Duration current_total_time() const;
    void print_time(Phase phase, const std::string& name) const;
    void print_stats() const;
    void write_stats() const;
private:
    const Options& options;
    std::vector<std::optional<TimePoint>> phase_latest_starttime;
    std::vector<Duration> durations;

public:

    static double cpuTime();
    static double childCpuTime();
    static double memUsage();
    static double peakMemUsage();

    //extra function to add time from child process used for clique to total times
    void add_clique_time(double time);

    // #########  all the variables and data points we want to observe and store during the algorithm  ####

    bool solved = false;

    //any lower and upper bounds + time to reach them
    int lower_bound;
    Duration tt_lower_bound;
    int clique_lower_bound;
    int mc_lower_bound;
    int upper_bound;
    Duration tt_upper_bound;
    int heuristic_bound;

    //preprocessing information
    std::tuple<int, int, double> original_graph_size;
    std::tuple<int, int, double> reduced_graph_size;
    int num_removed_small_degree = 0;
    int num_removed_dominated = 0;
    bool solved_in_preprocessing = false;

    //more fine-grained stats about the size of the instance
    int num_vars = 0;
    int num_clauses = 0;
    int num_sij_vars = 0;
    int num_cj_vars = 0;
    int num_transitivity_clauses = 0;
    int num_vars_at_most_k = 0;
    int num_clauses_at_most_k = 0;
    int num_removable_cj = 0;
    //at this point we also write peak memory usage!

    //cegar stats
    int num_cegar_iterations = 0;
    int num_total_conflicts = 0;
    std::vector<int> store_num_conflicts;

    //struct to store data for a bound
    using BoundInfo = std::tuple<
        //bound, satisfiable, total_time, cegar_time, solver_time, cegar_iterations, cegar_conflicts,
          int,   bool,        double,     double,     double,      int,              int,
        //max_level, assignments, decisions, backtracks, propagations, reason clauses, external clauses
          int,       long,        long,      long,       long,         long,           long
    >;
    BoundInfo produce_bound_info(int b, bool sat) {
        return {b, sat, current_total_time().count(), duration_of(CEGAR).count(),
                duration_of(SatSolver).count(), num_cegar_iterations, num_total_conflicts,
                prop_max_level, prop_num_assignments, prop_num_decisions, prop_num_backtracks,
                prop_num_propagations, prop_num_reason_clauses, prop_num_external_clauses};
    };

    //maps a bound to the number of iterations/conflicts/time it required
    std::vector<BoundInfo> bound_information;

    //Cadical propagator general statistics, stuff like max level, decisions, propagations...
    int prop_max_level = 0;
    long long prop_num_assignments = 0;
    long long prop_num_decisions = 0;
	std::vector<int> prop_node_depth_history; //tracks how often a node of level i was visited
    long long prop_num_backtracks = 0;
    std::vector<int> prop_backtrack_size; //stores size of jumps in backtracks
    //interesting statistic but vector becomes too large, only use this for testing on smaller instances
    std::vector<int> prop_detailed_backtrack_list; //stores old_lvl, new_lvl, bt_reason as adjacent entries
    int backtrack_reson = 0; //0 for cadical (likely conflict), 1 for clique and 2 for mycielsky bound
    long long prop_num_propagations = 0;
    long long prop_num_reason_clauses = 0;
    long long prop_num_external_clauses = 0;

    long long prop_num_clique_computations = 0;
    long long prop_num_maximal_cliques_computed = 0;
    long long prop_num_tight_cliques_computed = 0;
    long long prop_num_clique_successes = 0;
	std::vector<int> prop_clique_pruning_level;  //tracks for level i how often pruning was successful
    //collect for num_colors - i how often mycielsky was called and how often it was successful
    std::vector<int> mycielsky_calls = std::vector<int>(1, 0);
    std::vector<int> mycielsky_sucesses = std::vector<int>(1, 0);
	std::vector<int> prop_myc_pruning_level;
    long long prop_num_dominated_vertex_decisions = 0;
    long long prop_positive_prunings = 0;
    std::vector<int> prop_positive_pruning_level;
    long long prop_negative_prunings = 0;
    std::vector<int> prop_negative_pruning_level;

    //stats for coloring heuristic success, only interested in how much time could theoretically have been saved
    bool heuristic_found_coloring = false;
    //when the improved coloring was found
    Duration heuristic_coloring_time_point;
    //how much time would have been saved if stopping when heuristic found coloring
    Duration heuristic_theoretical_time_improvement = Duration(0);

};
std::ostream &operator<<(std::ostream &s, Duration duration);
template<typename T>
std::ostream &operator<<(std::ostream &s, const std::vector<T> &vec);
template<typename T>
std::string vec2str(const std::vector<T> &vec);
std::vector<int> truncate(std::vector<int> vec); //remove trailing zeros
//helper function to print tuples
template <size_t n, typename... T>
std::enable_if_t<(n >= sizeof...(T))> print_tuple(std::ostream&, const std::tuple<T...>&){} //end of recusion, do nothing
template <size_t n, typename... T>
std::enable_if_t<(n < sizeof...(T))> print_tuple(std::ostream& os, const std::tuple<T...>& tup); //print current index and call to print next
template <typename... T>
std::ostream& operator<<(std::ostream& os, const std::tuple<T...>& tup); //put tuple print into ostream


#endif //INCSATGC_STATISTICS_H
