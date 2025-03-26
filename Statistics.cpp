#include "Statistics.h"


Statistics::Statistics(const Options &options)
    : options(options)
{
    phase_latest_starttime.resize((static_cast<int>(Phase::Total)) + 1, {});  // Total is the last enum value
    durations.resize((static_cast<int>(Phase::Total)) + 1, {});
    start_phase(Phase::Total);
}

void Statistics::start_phase(const Statistics::Phase phase) {
    //set latest time point for started phase to late get the duration in end_phase
    if (!phase_latest_starttime.at(static_cast<int>(phase)).has_value()) {
        phase_latest_starttime.at(static_cast<int>(phase)) = cpuTime();
    }
}

void Statistics::end_phase(const Statistics::Phase phase) {
    //if phase has a latest start-time, use it to compute duration and end phase by resetting latest start-time
    auto start = phase_latest_starttime.at(static_cast<int>(phase));
    if (not start.has_value()) {
        //phase was not active and thus nothing done
        return;
    }
    //else
    //phase was active and is ended
    Duration duration = Duration(cpuTime() - start.value());
    durations.at(static_cast<int>(phase)) += duration;
    phase_latest_starttime.at(static_cast<int>(phase)).reset();

    //when ending a phase, make sure that previous phases are also stopped, to ensure good timekeeping when program is interrupted
    switch (phase) {
        case Total:
            end_phase(Preprocessing);
            end_phase(Algorithm);
            break;
        case Preprocessing:
            end_phase(PreprocessingClique);
            end_phase(PreprocessingMycielsky);
            end_phase(PreprocessingReductions);
            end_phase(PreprocessingInitialColoring);
            break;
        case Algorithm:
            end_phase(Preprocessing);
            end_phase(BuildEncoding);
            end_phase(SatSolver);
            break;
        case SatSolver:
            end_phase(PropagatorAssignment);
            end_phase(PropagatorBacktrack);
            end_phase(PropagatorDecide);
            end_phase(PropagatorComputeCliques);
            end_phase(PropagatorCliqueClauses);
            end_phase(PropagatorMycielskyClauses);
            end_phase(PropagatorPositivePruning);
            end_phase(PropagatorNegativePruning);
            end_phase(TestColorTime);
            break;
        default:
            break;
    }
}

Duration Statistics::duration_of(const Statistics::Phase phase) const {
    return durations.at(static_cast<int>(phase));
}

double Statistics::duration_in_sec(const Statistics::Phase phase) const {
    //return double of number of seconds up to millisecond precision
    return std::round(durations.at(static_cast<int>(phase)).count() * 1000) / 1000;
}

Duration Statistics::current_total_time() const {
    auto start = phase_latest_starttime.at(static_cast<int>(Total));
    assert(start.has_value());
    if(start.has_value()) {
        return Duration(cpuTime() - start.value()) + Duration(childCpuTime()); //also count time from children
    }
    return Duration{};
}

void Statistics::print_time(const Phase phase, const std::string& name) const {
    std::cout << "c Stats: Time " << name << ": " << duration_of(phase) << "\n";
}

void Statistics::print_stats() const {
    std::cout << "c #################################\n";
    std::cout << "c Stats: Instance " << options.filename << "\n";
    std::cout << "c Stats: Memory " << peakMemUsage() << " mb\n";
    //nice format printing of times for all phases
    print_time(Total,                       "┌ Total");
    print_time(Preprocessing,               "│ ┌ Preprocessing");
    print_time(PreprocessingClique,         "│ │ - Clique");
    print_time(PreprocessingMycielsky,      "│ │ - Mycielsky");
    print_time(PreprocessingReductions,     "│ │ - Reductions");
    print_time(PreprocessingInitialColoring,"│ └ - Initial coloring");
    print_time(Algorithm,                   "│ ┌ Algorithm");
    print_time(BuildEncoding,               "│ │ - Building encoding");
    print_time(BuildAtMostK,                "│ │ - Building at most k constraints");
    print_time(CEGAR,                       "│ │ - CEGAR");
    // if (options.encoding != Options::ZykovPropagator) {
    print_time(SatSolver,                   "└ └ - Sat solving");
    // }
    // else {
    // print_time(SatSolver,                   "│ │ ┌ Sat solving");
    // print_time(PropagatorAssignment,        "│ │ │ - Assign");
    // print_time(PropagatorBacktrack,         "│ │ │ - Backtrack");
    // print_time(PropagatorDecide,            "│ │ │ - Decide");
    // print_time(PropagatorComputeCliques,    "│ │ │ - Compute Cliques");
    // print_time(PropagatorCliqueClauses,     "│ │ │ - Clique Clauses");
    // print_time(PropagatorMycielskyClauses,  "│ │ │ - Mycielsky Clauses");
    // print_time(PropagatorPositivePruning,   "│ │ │ - Positive Pruning");
    // print_time(PropagatorNegativePruning,   "└ └ └ - Negative Pruning");
    // if (options.zykov_coloring_algorithm != Options::None) {
    // print_time(TestColorTime,               "(testing) - Coloring time");
    // }
    // }

    if (num_cegar_iterations > 0){
        std::cout << "c Stats: Cegar " << num_total_conflicts << " conflicts in "
                  << num_cegar_iterations << " iterations, taking overall " << duration_of(CEGAR) << "\n";
    }

    if(options.strategy == Options::SingleK){
        assert(options.specific_num_colors.has_value());
        std::cout << "c Stats: Instance is " << (upper_bound <= options.specific_num_colors ? "SAT" :
                                                     (lower_bound > options.specific_num_colors ? "UNSAT" : "UNKNOWN"))
                  << " for K = " << options.specific_num_colors.value() << "\n";
    }
    else{
        if(lower_bound == upper_bound){
            std::cout << "c Stats: Instance has chromatic number "<< lower_bound << "\n";
        }
        else{
            std::cout << "c Stats: Instance has lower bound "<< lower_bound << " and upper bound " << upper_bound << "\n";
        }
    }
    if(options.encoding != Options::ZykovPropagator) {
        return;
    }
    if(options.verbosity <= Options::Normal) { //only print these stats as extra information
        return;
    }
    double total_sat_time = duration_in_sec(SatSolver);
    std::cout << "Propagator stats: \nmax level "  << prop_max_level
              << "\nassignments "  << prop_num_assignments << " (" << prop_num_assignments/total_sat_time << " /sec)"
              << "\ndecisions "  << prop_num_decisions << " (" << prop_num_decisions/total_sat_time << " /sec)" //<< " history " << truncate(prop_node_depth_history)
              << "\nbacktracks "  << prop_num_backtracks << " (" << prop_num_backtracks/total_sat_time << " /sec)" //<< " history " << truncate(prop_backtrack_size) //<< " history " << prop_detailed_backtrack_list
              << "\npropagations "  <<  prop_num_propagations << " (" << prop_num_propagations/total_sat_time << " /sec)"
              << "\nreason clauses "  << prop_num_reason_clauses << "\nexternal clauses "  << prop_num_external_clauses
              << "\nnum clique computations " << prop_num_clique_computations << "\nnum maxcliques computed " << prop_num_maximal_cliques_computed
              << "\nnum tight cliques computed " << prop_num_tight_cliques_computed
              << "\nclique successes " << prop_num_clique_successes //<< " history " << truncate(prop_clique_pruning_level)
              << "\nmycielsky calls " << mycielsky_calls
              << "\nmycielsky successes " << mycielsky_sucesses //<< " history " << truncate(prop_clique_pruning_level)
              << "\ndominated vertices " << prop_num_dominated_vertex_decisions
              << "\npositive prunings " << prop_positive_prunings //<< " history " << truncate(prop_positive_pruning_level)
              << "\nnegative prunings " << prop_negative_prunings //<< " history " << truncate(prop_negative_pruning_level)
              << "\n";
              if (options.zykov_coloring_algorithm != Options::None) {
              std::cout << "coloring heuristic time saved " << heuristic_theoretical_time_improvement << "\n";
              }
}

void Statistics::write_stats() const {
    //first row of csv file naming all the data entries. Take care to keep the names and entries synchronised
    std::stringstream csv_header;
    csv_header << "instance;" "full cmd;" "solved;"
               //time stats
               "total time;" "preprocessing time;" "clique time;" "mycielsky time;" "reductions time;" "initial coloring time;"
               "algorithm time;" "encoding time;" "at most k time;" "cegar time;" "solver time;"
               //algorithm information stats
               "lower bound;" "tt_lb;" "clique bound;" "mycielsky bound;"
               "upper bound;" "tt_ub;" "heuristic bound;"
               "original n;" "original m;" "original d;" "reduced n;" "reduced m;" "reduced d;"
               "removed smalldegree;" "removed dominated;" "solved in preprocessing;"
               "variables;" "clauses;" "s_ij variables;" "cj variables;" "transitivity clauses;"
               "at most k variables;" "at most k clauses;" "removed cj;" "peak memory usage;"
               "cegar iterations;" "conflicts;" "conflicts per iteration;"
               "bound information history;"
               //propagator stats
               "max level;" "assignments;" "decisions;" "decision levels;" "backtracks;" "backtrack size;" //"backtrack list;"
               "propagations;" "reason clauses;" "external clauses;"
               "cliques computations;" "maximal cliques computed;" "tight cliques computed;"
               "cliques successes;" "clique levels;"
               "mycielsky calls;" "mycielsky successes;" "mycielsky levels;"
               "dominated vertex decisions;" "positive prunings;" "positive pruning levels;"
               "negative prunings;" "negative pruning levels"
               <<(options.zykov_coloring_algorithm != Options::None ? ";heuristic time improvement" : "")  //optional stats that are not always reported
               <<(options.enable_detailed_backtracking_stats ? ";detailed backtrack stats" : "")<<
               "\n";
    // If the csv file does not exist, create it and write first row, otherwise open and append
    std::ofstream csv_file;
    if (not std::ifstream(options.stats_csvfile)) {
        csv_file.open(options.stats_csvfile);
        csv_file << csv_header.str();
    } else {
        csv_file.open(options.stats_csvfile, std::ios_base::app);
    }
    //write the actual data points in the next row of the csv file
    csv_file << options.filename << ";" << options.full_cmd << ";" << solved << ";"
            //time stats
            << duration_in_sec(Total) << ";" << duration_in_sec(Preprocessing) << ";"
            << duration_in_sec(PreprocessingClique) << ";" << duration_in_sec(PreprocessingMycielsky) << ";"
            << duration_in_sec(PreprocessingReductions) << ";" << duration_in_sec(PreprocessingInitialColoring) << ";"
            << duration_in_sec(Algorithm) << ";" << duration_in_sec(BuildEncoding) << ";"
            << duration_in_sec(BuildAtMostK) << ";" << duration_in_sec(CEGAR) << ";" << duration_in_sec(SatSolver) << ";"
            //algorithm information stats
            << lower_bound << ";" << tt_lower_bound.count() << ";" << clique_lower_bound << ";" << mc_lower_bound << ";"
            << upper_bound << ";" << tt_upper_bound.count() << ";" << heuristic_bound << ";"
            << std::get<0>(original_graph_size) << ";" << std::get<1>(original_graph_size) << ";" << std::get<2>(original_graph_size) << ";"
            << std::get<0>(reduced_graph_size) << ";" << std::get<1>(reduced_graph_size) << ";" << std::get<2>(reduced_graph_size) << ";"
            << num_removed_small_degree << ";" << num_removed_dominated << ";" << solved_in_preprocessing << ";"
            << num_vars << ";" << num_clauses << ";" << num_sij_vars << ";" << num_cj_vars << ";" << num_transitivity_clauses << ";"
            << num_vars_at_most_k << ";" << num_clauses_at_most_k << ";" << num_removable_cj << ";" << peakMemUsage() << ";"
            << num_cegar_iterations << ";" << num_total_conflicts << ";" << store_num_conflicts << ";"
            << bound_information << ";"
            //propagator stats
            << prop_max_level << ";" << prop_num_assignments << ";"
            << prop_num_decisions << ";" << truncate(prop_node_depth_history) << ";"
            << prop_num_backtracks  << ";" << truncate(prop_backtrack_size) << ";"
            << prop_num_propagations << ";" << prop_num_reason_clauses << ";" << prop_num_external_clauses << ";"
            << prop_num_clique_computations << ";" << prop_num_maximal_cliques_computed << ";" << prop_num_tight_cliques_computed << ";"
            << prop_num_clique_successes << ";" << truncate(prop_clique_pruning_level) << ";"
            << mycielsky_calls << ";" << mycielsky_sucesses << ";" << truncate(prop_myc_pruning_level) << ";"
            << prop_num_dominated_vertex_decisions << ";" << prop_positive_prunings << ";" << truncate(prop_positive_pruning_level) << ";"
            << prop_negative_prunings << ";" << truncate(prop_negative_pruning_level)
            << (options.zykov_coloring_algorithm != Options::None ? ";"+std::to_string(heuristic_theoretical_time_improvement.count()) : "")
            << (options.enable_detailed_backtracking_stats? ";"+vec2str(prop_detailed_backtrack_list) : "")
            << "\n";
}


double Statistics::cpuTime() {
#ifdef __unix__
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return static_cast<double>(ru.ru_utime.tv_sec)
           + static_cast<double>(ru.ru_utime.tv_usec) / 1000000;
#else
    return static_cast<double>(std::clock()) /  CLOCKS_PER_SEC;
#endif
}

double Statistics::childCpuTime() {
#ifdef __unix__
    struct rusage ru;
    getrusage(RUSAGE_CHILDREN, &ru);
    return static_cast<double>(ru.ru_utime.tv_sec)
           + static_cast<double>(ru.ru_utime.tv_usec) / 1000000;
#else
    return 0.0;
#endif
}

double Statistics::memUsage() {
    double memory_usage = 0.0;
#ifdef __unix__
    std::string line;
    std::ifstream stat_stream("/proc/self/stat", std::ios_base::in);
    if (stat_stream) {
        std::string dummy;
        long rss; //resident set size
        stat_stream >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy
            >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy
            >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy
            >> dummy >> dummy >> rss; //only want rss
        stat_stream.close();

        long page_size_kb = sysconf(_SC_PAGE_SIZE)
            / 1024; // in case x86-64 is configured to use 2MB pages
        memory_usage = rss * page_size_kb;
        memory_usage = memory_usage / 1024;//want mb
        memory_usage = std::round(memory_usage * 1000) / 1000; //round to two places
    }
#endif
    return memory_usage;
}

double Statistics::peakMemUsage() {
#ifdef __unix__
    struct rusage u;
    getrusage (RUSAGE_SELF, &u);
    //max resident set size in kb, divide by 1024 to get mb
    return static_cast<double>(u.ru_maxrss) / 1024;
#else
    return 0.0;
#endif
}


void Statistics::add_clique_time(const double time){
    Duration duration(time);
    //have to increase Total, Preprocessing and PreprocessingClique
    durations.at(static_cast<int>(Total)) += duration;
    durations.at(static_cast<int>(Preprocessing)) += duration;
    durations.at(static_cast<int>(PreprocessingClique)) += duration;
}



std::ostream &operator<<(std::ostream &s, Duration duration) {
    auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
    duration -= hours;
    if (hours.count()) {
        s << hours.count() << "h ";
    }

    auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration);
    duration -= minutes;
    if (minutes.count()) {
        s << minutes.count() << "m ";
    }

    auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration);
    duration -= seconds;
    s << seconds.count() << "s ";

    auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(duration);
    duration -= milliseconds;
    s << milliseconds.count() << "ms";
    //stop at milliseconds, we don't need any more precision than that
    return s;
}


template<typename T>
std::ostream &operator<<(std::ostream &s, const std::vector<T> &vec) {
    s << "[";
    std::string sep;
    for(T el : vec){
        s << sep << el;
        sep = ", ";
    }
    return s << "]";
}

template<typename T>
std::string vec2str(const std::vector<T> &vec) {
    std::stringstream s;
    s << "[";
    std::string sep;
    for(T el : vec){
        s << sep << el;
        sep = ", ";
    }
    s << "]";
    return s.str();
}


std::vector<int> truncate(std::vector<int> vec) {
    const auto it = std::find_if(vec.rbegin(), vec.rend(), [](const int x) { return x != 0; });
    vec.erase(it.base(), vec.end()); //return vector but keep last 0
    return vec;
}



template<size_t n, typename ... T>
std::enable_if_t<(n < sizeof...(T))> print_tuple(std::ostream &os, const std::tuple<T...> &tup) {
    if (n != 0) {
        os << ", ";
    }
    os << std::get<n>(tup);
    print_tuple<n+1>(os, tup);
}

template<typename ... T>
std::ostream & operator<<(std::ostream &os, const std::tuple<T...> &tup) {
    os << "[";
    print_tuple<0>(os, tup);
    return os << "]";
}
