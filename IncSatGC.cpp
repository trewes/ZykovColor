#include "IncSatGC.h"


UpperTriangle::UpperTriangle() : dimension(0){
    //default constructor
}

UpperTriangle::UpperTriangle(int dimension, int value) : dimension(dimension){
    //number of element of upper triangular matrix without diagonal
    int num_elements = (dimension-1) * dimension / 2;
    sij_indices.resize(num_elements, value);
    index_to_ij.resize(num_elements);
}


IncSatGC::IncSatGC(const char *filename, Options options) :
        IncSatGC(Graph::Graph(filename), std::move(options)) //build graph and call other constructor with that
{}

IncSatGC::IncSatGC(Graph::Graph in_graph, Options options_) :
        graph(std::move(in_graph)),
        num_vertices(static_cast<int>(graph.ncount())),
        options(std::move(options_)),
        stats(Statistics(options)),
        lower_bound(1 + (graph.ecount() > 0)),
        clique_lower_bound(1 + (graph.ecount() > 0)),
        mc_lower_bound(1),
        frac_lower_bound(0),
        upper_bound(num_vertices),
        heuristic_bound(num_vertices),
        complement_graph_adjacency({}),
        sij_indices({}),
        c_indices({}),
        encoder(),
        encoder_assumptions({}),
        cj_literals({}),
        num_removable_cj(0),
        is_removable_cj({}),
        current_best_coloring({}),
        solved_in_preprocessing(false)
{
    if(options.verbosity >= Options::Normal) {
        options.print_header();
        options.print();
    }
    //just to keep stats synchronised in extreme edge cases
    stats.lower_bound = lower_bound;
    stats.clique_lower_bound = clique_lower_bound;
    stats.mc_lower_bound = mc_lower_bound;
    stats.frac_lower_bound = frac_lower_bound;
    stats.upper_bound = upper_bound;
    stats.heuristic_bound = heuristic_bound;
    new_SAT_solver();

    //for writing cleanup on exit/interrupt
    INSTANCE = this;
    register_write_cleanup_on_exit();
}


int IncSatGC::run() {
    stats.start_phase(Statistics::Total);

    if(not options.disable_preprocessing){
        preprocessing();
        if (solved_in_preprocessing){
            assert(lower_bound == upper_bound);
            std::cout << "Result: Chromatic number of " << lower_bound << " was determined in pre-processing.\n";
            stats.end_phase(Statistics::Total);
            write_and_cleanup();
            return lower_bound;
        }
    }

    //main execution of algorithm to find chromatic number
    int chromatic_number = compute_chromatic_number();

    stats.end_phase(Statistics::Total);
    write_and_cleanup();
    return chromatic_number;
}



void IncSatGC::preprocessing(){
    stats.start_phase(Statistics::Preprocessing);
    stats.original_graph_size = {graph.ncount(), graph.ecount(), graph.density()};

    //get lower bounds on chromatic number and compute first initial coloring
    preprocessing_clique_bound();
    preprocessing_mycielsky_bound();
    preprocessing_initial_coloring();
    //check if that already solved the graph
    if(lower_bound == upper_bound){
        solved_in_preprocessing = true;
        stats.solved = true;
        stats.solved_in_preprocessing = solved_in_preprocessing;
        return;
    }
    //use bound and dominated vertices to reduce graph
    if( options.reduce_graph ){
        preprocessing_reductions();
    }
    //find an initial coloring solution and upper bound
    if(not solved_in_preprocessing) {
        preprocessing_initial_coloring();
    }
    if(lower_bound == upper_bound){
        solved_in_preprocessing = true;
        stats.solved = true;
        stats.solved_in_preprocessing = solved_in_preprocessing;
        return;
    }
    //permute the graph if option is set
    if ((not solved_in_preprocessing) and options.use_clique_in_ordering) {
        preprocessing_clique_ordering();
    }
    //compute fractional chromatic number as initial bound (with a timeout)
    if (not solved_in_preprocessing) {
        frac_lower_bound = external_get_fractional();
        notify_lower_bound(frac_lower_bound);
        if (options.verbosity >= Options::Debug) {
            std::cout << "c Initial fractional lower bound of " << std::setprecision(25) << frac_lower_bound << std::setprecision(5)
                      << " in " << stats.duration_of(Statistics::PreprocessingFractional) <<  "\n";
        }
        if(lower_bound == upper_bound){
            solved_in_preprocessing = true;
            stats.solved = true;
            stats.solved_in_preprocessing = solved_in_preprocessing;
        }
    }

    if (options.verbosity >= Options::Normal) {
        std::cout << "c Found lb " << lower_bound << " and ub " << upper_bound <<
                  " at " << stats.current_total_time() << " and " << Statistics::memUsage() << " mb memory" ;
        if (options.verbosity >= Options::Verbose) {
            std::cout << " (clique lb " << clique_lower_bound << ", mycielsky lb " << mc_lower_bound
                      << ", fractional lb "<< frac_lower_bound << ", heuristic ub " << heuristic_bound << ")";
        }
        std::cout << ".\n";
    }

    stats.end_phase(Statistics::Preprocessing);
}

int IncSatGC::compute_chromatic_number() {
    stats.start_phase(Statistics::Algorithm);
    if(options.original_paper_configuration) {
        original_paper_configuration();
        stats.end_phase(Statistics::Algorithm);
        assert(lower_bound == upper_bound);
        return upper_bound;
    }
    switch (options.encoding) {
        case Options::AssignmentEncoding:
            do_assignment_encoding();
            break;
        case Options::FullEncoding:
            do_full_encoding();
            break;
        case Options::FullMaxSAT:
            write_full_maxsat_encoding();
            break;
        case Options::CEGAR:
            do_cegar_approach();
            break;
        case Options::PartialOrderEncoding:
            do_partial_order_encoding();
            break;
        case Options::ZykovPropagator:
            do_zykov_propagator();
            break;
        default:
            throw std::runtime_error("Invalid option for encoding.");
    }
    stats.end_phase(Statistics::Algorithm);
    if(options.strategy == Options::SingleK) { //only SAT or UNSAT answer for chosen k
        assert(lower_bound == (options.specific_num_colors.value() + 1) or upper_bound == options.specific_num_colors);
        if(upper_bound <= options.specific_num_colors) {
            return 1; //SAT
        }
        if(lower_bound > options.specific_num_colors) {
            return 0; //UNSAT
        }
    }
    assert(lower_bound == upper_bound);
    return upper_bound;
}


void IncSatGC::new_SAT_solver() {
    if (solver != nullptr) {
        throw std::runtime_error("Solver should be NULL when initialising new solver.");
    }
    switch (options.solver) {
        case Options::CaDiCaL:
            solver = std::make_shared<CaDiCaLAdaptor::Solver>();
            break;
        case Options::Glucose:
            solver = std::make_shared<NSPACE::Solver>();//using SimpSolver breaks things but I don't know why
            solver->verbosity = 0;
//            solver->setIncrementalMode();//setIncremental breaks things, incremental solve works anyway
            break;
        default:
            throw std::runtime_error("Invalid option for solver.");
    }
    if (solver == nullptr) {
        throw std::runtime_error("Unable to create solver.");
    }
}

void IncSatGC::reset_SAT_solver() {
    if (solver == nullptr) {
        throw std::runtime_error("Solver should not be NULL when resetting solver.");
    }
    solver = nullptr;
    //propagators are reset automatically when a new one is assigned since unqie_ptr is used
    new_SAT_solver();
}

void IncSatGC::add_vars(const int num_vars) {
    //we have to add vars to base glucose solver in any case, since that count is used in the encoding

    if (options.solver == Options::CaDiCaL){
        std::shared_ptr<CaDiCaLAdaptor::Solver> cast_solver = std::dynamic_pointer_cast<CaDiCaLAdaptor::Solver>(solver);
        if (not cast_solver) { //solver ptr was not of type CaDiCaLAdaptor::Solver
            throw std::runtime_error("Solver* was of the wrong type.");
        }
        //added functionality to add multiple variables at once in ExtendSolvers
        cast_solver->newVars(num_vars);
    }
    else if (options.solver == Options::Glucose) {
        //for base glucose solver, add vars one by one
        for (int i = 0; i < num_vars; ++i) {
            solver->newVar();
        }
    }
    else {
        throw std::runtime_error("Invalid option for solver.");
    }
}


int IncSatGC::get_num_vars() const {
    if (options.solver == Options::CaDiCaL) {
        std::shared_ptr<CaDiCaLAdaptor::Solver> cast_solver = std::dynamic_pointer_cast<CaDiCaLAdaptor::Solver>(solver);
        if (not cast_solver) { //solver ptr was not of type CaDiCaLAdaptor::Solver
            throw std::runtime_error("Solver* was of the wrong type.");
        }
        return cast_solver->solver.vars();
    }
    if (options.solver == Options::Glucose) {
        return solver->nVars();
    }
    throw std::runtime_error("Invalid option for solver.");
}

int IncSatGC::get_num_clauses() const {
    if (options.solver == Options::CaDiCaL) {
        std::shared_ptr<CaDiCaLAdaptor::Solver> cast_solver = std::dynamic_pointer_cast<CaDiCaLAdaptor::Solver>(solver);
        if (not cast_solver) { //solver ptr was not of type CaDiCaLAdaptor::Solver
            throw std::runtime_error("Solver* was of the wrong type.");
        }
        return cast_solver->num_clauses;
    }
    if (options.solver == Options::Glucose) {
        return solver->nClauses();
    }
        throw std::runtime_error("Invalid option for solver.");
}

bool IncSatGC::run_solver() {
    stats.start_phase(Statistics::SatSolver);

    if (options.solver == Options::CaDiCaL) {
        std::shared_ptr<CaDiCaLAdaptor::Solver> cast_solver = std::dynamic_pointer_cast<CaDiCaLAdaptor::Solver>(solver);
        if (not cast_solver) { //solver ptr was not of type CaDiCaLAdaptor::Solver
            throw std::runtime_error("Solver* was of the wrong type.");
        }

        cast_solver->assume(encoder_assumptions);
        int result = cast_solver->solver.solve();
        if (result == CaDiCaL::Status::UNKNOWN) {//solver inconclusive, cadical returns code 0
            throw std::runtime_error("Problem unsolved by CaDiCal.");
        }
        stats.end_phase(Statistics::SatSolver);
        return result == CaDiCaL::Status::SATISFIABLE;//cadical code for satisfiable
    }
    if (options.solver == Options::Glucose) {
        bool result = solver->solve(encoder_assumptions);
        stats.end_phase(Statistics::SatSolver);
        return result;
    }
        throw std::runtime_error("Invalid option for solver.");
}

void IncSatGC::write_cnf() {
    assert(options.strategy == Options::SingleK and options.specific_num_colors.has_value());

    std::string cnf_name = options.filename;
    switch (options.encoding) {
        case Options::AssignmentEncoding:
            cnf_name += "_assignment_encoding";
            break;
        case Options::PartialOrderEncoding:
            cnf_name += "_partial_order_encoding";
            break;
        case Options::FullEncoding:
            cnf_name += "_full_zykov_encoding";
            break;
        case Options::FullMaxSAT:
        case Options::CEGAR:
        case Options::ZykovPropagator:
        default:
            std::cout << "Writing of cnf not supported for the given encoding\n";
            return;
    }
    cnf_name += "_"+std::to_string(options.specific_num_colors.value());
    cnf_name += ".cnf";


    if (options.solver == Options::CaDiCaL) {
        std::shared_ptr<CaDiCaLAdaptor::Solver> cast_solver = std::dynamic_pointer_cast<CaDiCaLAdaptor::Solver>(solver);
        if (not cast_solver) { //solver ptr was not of type CaDiCaLAdaptor::Solver
            throw std::runtime_error("Solver* was of the wrong type.");
        }
        cast_solver->solver.write_dimacs(cnf_name.c_str());
    }
    else if (options.solver == Options::Glucose){
        solver->toDimacs(cnf_name.c_str());
    }
    else {
        throw std::runtime_error("Invalid option for solver.");
    }
}

Model IncSatGC::get_model() const {
    Model model = Model(get_num_vars());
    if (options.solver == Options::CaDiCaL) {
        std::shared_ptr<CaDiCaLAdaptor::Solver> cast_solver = std::dynamic_pointer_cast<CaDiCaLAdaptor::Solver>(solver);
        if (not cast_solver) { //solver ptr was not of type CaDiCaLAdaptor::Solver
            throw std::runtime_error("Solver* was of the wrong type.");
        }
        for (int i = 1; i <= cast_solver->solver.vars(); ++i) {
            if(cast_solver->solver.val(i) > 0){
                model.set(i - 1); //cadical indices start at 1
            }
        }
    }
    else if (options.solver == Options::Glucose){
        for (int i = 0; i < solver->nVars(); ++i) {
            if(solver->modelValue(i) == l_True){
                model.set(i);
            }
        }
    }
    else {
        throw std::runtime_error("Invalid option for solver.");
    }
    return model;
}

bool IncSatGC::get_model_ij(const Model &model, int i, int j) const {
    //computes (sij_indices.get(i,j) == -1) ? false : model[sij_indices.get(i,j)];
    //i.e. if there is an edge, return false, otherwise return what the solution assigned
    //can logically be simplified to the following
    return (sij_indices.get(i, j) != -1) && model[sij_indices.get(i, j)];
}


std::vector<Graph::VertexType> IncSatGC::external_get_clique(bool write_graph_get_clique) {
    std::vector<Graph::VertexType> CliSAT_clique;
    namespace bp = boost::process::v1;

    //first check that external CliSAT binary for cliques exists
    std::filesystem::path CliSAT_path;
#ifdef CLISAT_BINARY_PATH
    CliSAT_path = CLISAT_BINARY_PATH;
    if (not std::filesystem::exists(CliSAT_path)) {
        throw std::runtime_error("The external CliSAT binary does not exist at the specified path.");
    }
#else
    std::cout << "WARNING: CliSAT binary was not specified so the algorithm does not use an initial clique, "
                 "likely degrading performance significantly.\n";
    return {};
#endif

    std::string tmp_name;
    if(write_graph_get_clique){
        //write current problem graph to file, so we can run CliSAT
        tmp_name = "tmp_reduced_" + std::string(options.filename);
        graph.write_dimacs(tmp_name);
    }

    std::vector<std::string> args(4);
    args[0] = (write_graph_get_clique ? tmp_name : options.filepath);
    args[1] = "1";//add time limit of 1 second
    args[2] = "1";//and choose method 1 for CliSAT
    args[3] = "1";//enable AMTS heuristic for CliSAT

    std::string cmd_output;
    bp::ipstream out_stream; // stream for command output

    auto start = Statistics::childCpuTime();
    // execute the command and redirect the output to out_stream
    bp::child c(CliSAT_path.string(), args, bp::std_out > out_stream);
    // wait for the process to finish
    c.wait();
    //add child process time to stats
    stats.add_clique_time(Statistics::childCpuTime() - start);

    // read the command's output and store all lines
    std::vector<std::string> command_output_lines;
    while(std::getline(out_stream, cmd_output)){
        command_output_lines.push_back(cmd_output);
    }
    std::string clique_output = *(command_output_lines.rbegin());
    clique_output = clique_output.substr(7,clique_output.length());
    //extract size of clique and truncate string
    int clique_size = 0;
    std::size_t found = clique_output.find('[');
    if (found != std::string::npos) {
        //new string of only number in brackets
        std::string out_size_str = clique_output.substr(found + 1, clique_output.find(']') - found - 1);
        clique_size = std::stoi(out_size_str);
        // Truncate original string
        clique_output = clique_output.substr(0, found);
    }
    CliSAT_clique.reserve(clique_size);
    std::istringstream iss(clique_output);
    Graph::VertexType v;
    while (iss >> v){
        CliSAT_clique.push_back(v); //their vertex numbering starts at 0
    }

    //check whether clique found was optimal
    bool clique_is_optimal = true;
    for (const std::string& line : command_output_lines) {
        if(line.find("exiting on time out") != std::string::npos){
            clique_is_optimal = false;
            break;
        }
    }
    //print results
    if(options.verbosity >= Options::Verbose) {
        std::cout << "c Clique: " << (write_graph_get_clique ? "[In reduced Graph] " : "") << "Found "
                  << (clique_is_optimal ? "an optimal" : "a non-optimal")
                  << " clique of size " << clique_size ;
        if(options.verbosity >= Options::Debug) {
            std::cout << ": ";
            for (Graph::VertexType clique_vertex: CliSAT_clique) {
                std::cout << clique_vertex << " ";
            }
        }
        std::cout << "\n";
    }
    assert(static_cast<int>(CliSAT_clique.size()) == clique_size);

    if(write_graph_get_clique){
        //delete written temporary graph
        std::filesystem::remove(tmp_name);
    }
    return CliSAT_clique;
}

int IncSatGC::mycielsky_extension_lb(SubGraph H) const {
    //function to extend an initial subgraph to larger generalised mycielsky graph
    // returned is the size of H plus the number of successful extensions

    //get graph representation as vector of bitsets
    std::vector<Bitset> matrix(num_vertices, Bitset(num_vertices));
    auto elist = graph.elist();
    for(int i = 0; i < graph.ecount(); i++){
        matrix[elist[2 * i]].set( elist[2 * i + 1]);
        matrix[elist[2 * i + 1]].set( elist[2 * i]);
    }
    //start algorithm with H being a clique, so initial k is the size of the clique
    int initial_k = static_cast<int>(H.size());
    int result_k = initial_k;
    Bitset W = Bitset(num_vertices);
    while (static_cast<int>(H.size()) < num_vertices){
        //initialize W=V and S_v = {v}
        W.set();
        std::map< Graph::VertexType, Bitset > vecS;
        for (const auto& pair : H) {
            Graph::VertexType v = pair.first;
            vecS[v] = Bitset(num_vertices);
            vecS[v].set(v);
        }
        //first loop computing S_v and intersecting W
        for (const auto& pair : H) {
            Graph::VertexType v = pair.first;
            for (Graph::VertexType u = 0; u < num_vertices; ++u) {
                bool contains = (matrix[u] & H[v]) == H[v];
                if (contains) {
                    vecS[v].set(u);
                }
            }
            Bitset NG_Sv = Bitset(num_vertices);
            for (auto u = vecS[v].find_first(); u != Bitset::npos; u = vecS[v].find_next(u)){
                NG_Sv |= matrix[u]; //union
            }
            W &= NG_Sv; //intersection
            if(W.none()){ //W already empty, can stop here
                break;
            }
        }
        //case that W is not empty, can extend subgraph and increase lower bound
        if(W.any()){
            result_k++;
            Graph::VertexType w = W.find_first();
            auto H_old = H;
            for (const auto& pair : H_old) {
                Graph::VertexType v = pair.first;
                Bitset  NG_w_intersect_S_v = Bitset(num_vertices);
                NG_w_intersect_S_v = matrix[w] & vecS[v];
                Graph::VertexType u = NG_w_intersect_S_v.find_first();
                //insert vertex u from the intersection and insert edges
                //the followings inserts the element with empty bitset if it doesn't exist, otherwise finds iterator
                //either case, uses that iterator to set element in bitset
                //in u, insert w
                H.insert({u, Bitset(num_vertices)}).first->second.set(w);
                //in w, insert u
                H.insert({w, Bitset(num_vertices)}).first->second.set(u);
                for(Graph::VertexType vertex = H_old[v].find_first(); vertex != Bitset::npos; vertex = H_old[v].find_next(vertex)){
                    H[u].set(vertex);
                    H[vertex].set(u);
                }
            }
        }
        else {
            break;
        }
    }
    if(options.verbosity >= Options::Verbose){
        std::cout << "c Mycielsky: increased bound from initial " << initial_k << " to " << result_k << "\n";
    }
    return result_k;
}

bool IncSatGC::reduced_graph() {
    if(options.verbosity >= Options::Verbose) {
        std::cout << "c Preprocessing: original graph has " << graph.ncount() << " vertices, "
                  << graph.ecount() << " edges and density " << graph.density() * 100 << "%\n";
    }
    //preprocessing loop until graph doesn't change any further
    int original_num_nodes = graph.ncount();
    int num_nodes_previous_iteration = 0;
    while(num_nodes_previous_iteration != graph.ncount()){
        num_nodes_previous_iteration = graph.ncount();

        //peel graph with lower bound
        int num_removed_small_degree = graph.peel_graph(lower_bound);
        stats.num_removed_small_degree += num_removed_small_degree;

        if(graph.ncount() <= lower_bound) {
            //can return with lower bound, this already solves the graph
            notify_upper_bound(lower_bound);
            break;
        }

        //look for and remove dominated vertices
        int num_removed_dominated = graph.remove_dominated_vertices();
        stats.num_removed_dominated += num_removed_dominated;

        if (graph.ncount() <= lower_bound) {
            //can return with lower bound, this already solves the graph
            notify_upper_bound(lower_bound);
            break;
        }
    }
    //done with removing vertices, report change
    if(lower_bound == upper_bound) {
        solved_in_preprocessing = true;
        stats.solved = true;
        stats.solved_in_preprocessing = solved_in_preprocessing;
        if(options.verbosity >= Options::Normal) {
            std::cout << "c Preprocessing: Graph was solved in preprocessing.\n";
        }
    }
    else {
        if(options.verbosity >= Options::Verbose) {
            std::cout << "c Preprocessing: reduced graph has " << graph.ncount() << " vertices, "
                      << graph.ecount() << " edges and density " << graph.density() * 100 << "%";
            if(original_num_nodes != graph.ncount()){
                std::cout << ". removed " << stats.num_removed_small_degree << " small degree and "
                          << stats.num_removed_dominated << " dominated vertices";
            }
            std::cout << "\n";
        }
    }
    //return true if reduction happened, false if the graph is unchanged
    return original_num_nodes != graph.ncount();
}

void IncSatGC::preprocessing_clique_ordering() {
    //careful, this invalidates fields 'removed_vertices', 'removed_edges' and 'remap'
    //have to reverse permutation before using those fields again, i.e. before vertex recovery
    Graph::Permutation inverse_ordering = Graph::perm_inverse(graph_vertex_ordering);
    graph = graph.perm_graph(inverse_ordering);
    clique = Graph::apply_permutation(clique, inverse_ordering);
    //also apply permutation to current best/heuristic coloring, needed if we prove it to be optimal
    assert(not current_best_coloring.empty());
    auto tmp_coloring = current_best_coloring;
    for (int i = 0; i < static_cast<int>(current_best_coloring.size()); ++i) {
        current_best_coloring[inverse_ordering[i]] = tmp_coloring[i];
    }
}

void IncSatGC::preprocessing_initial_coloring() {
    stats.start_phase(Statistics::PreprocessingInitialColoring);
    graph_vertex_ordering.clear();
    heuristic_coloring = graph.dsatur(graph_vertex_ordering, clique);
    // graph_vertex_ordering = graph.max_connected_degree_ordering(clique);
    //coloring uses fewer colors than lower bound because of graph reduction, set upper bound = lower bound
    if(static_cast<int>(heuristic_coloring.size()) < lower_bound){
        assert(has_removed_vertices_in_reduction);
        assert(static_cast<int>(heuristic_coloring.size()) == lower_bound - 1);
        // add empty color class for recovered vertex later on
        heuristic_coloring.emplace_back();
    }
    notify_heuristic_ub(static_cast<int>(heuristic_coloring.size()));

    if (static_cast<int>(heuristic_coloring.size()) < heuristic_bound) {
        //transform and store heuristic coloring in current_best_coloring
        current_best_coloring.resize(num_vertices, NoColor);
        for (int cclass = 0; cclass < static_cast<int>(heuristic_coloring.size()); cclass++) {
            for (auto v : heuristic_coloring[cclass]) {
                current_best_coloring[v] = cclass;
            }
        }
        assert(std::count(current_best_coloring.begin(), current_best_coloring.end(), NoColor) == 0);
    }
    if(options.verbosity >= Options::Verbose){
        std::cout << "c Dsatur heuristic: computed upper bound of " << heuristic_coloring.size() << "\n";
    }
    stats.end_phase(Statistics::PreprocessingInitialColoring);
}

void IncSatGC::preprocessing_reductions() {
    stats.start_phase(Statistics::PreprocessingReductions);
    has_removed_vertices_in_reduction = reduced_graph();
    if(has_removed_vertices_in_reduction) {
        //update graph related numbers
        num_vertices = static_cast<int>(graph.ncount());
        if (not solved_in_preprocessing) {
            notify_upper_bound(num_vertices);
            heuristic_bound = num_vertices;
            //graph was reduced and some vertices removed, clique no longer valid in new graph
            // obtain clique from writing graph to file and running CliSAT again
            preprocessing_clique_bound(true);
            //also try to improve mycielsky bound again
            preprocessing_mycielsky_bound();
        }
        else {
            assert(lower_bound == upper_bound);
            //completely solved in reductions preprocessing, build full coloring later when recovering vertices
            heuristic_coloring = Graph::Coloring(lower_bound);
            if (graph.ncount()) {
                //there were still vertices left after reduction but less than lower_bound,
                // still have to assign them a color before recovering reduction
                assert(graph.ncount() <= lower_bound);
                for (int i = 0; i < graph.ncount(); ++i) {
                    heuristic_coloring[i].insert(i);
                }
            }
        }
    }
    stats.reduced_graph_size = {graph.ncount(), graph.ecount(), graph.density()};
    stats.end_phase(Statistics::PreprocessingReductions);
}

void IncSatGC::preprocessing_clique_bound(bool write_new_graph_and_get_clique) {
    stats.start_phase(Statistics::PreprocessingClique);
    clique = external_get_clique(write_new_graph_and_get_clique);
    notify_clique_lb(static_cast<int>(clique.size()));
    stats.end_phase(Statistics::PreprocessingClique);
}

void IncSatGC::preprocessing_mycielsky_bound() {
    if(options.use_mycielsky_lb) {
        stats.start_phase(Statistics::PreprocessingMycielsky);
        int mycielsky_bound = mycielsky_extension_lb(SubGraph({{0, Bitset(num_vertices)}})); //trivial subgraph for now
        notify_mycielsky_lb(mycielsky_bound);
        stats.end_phase(Statistics::PreprocessingMycielsky);
    }
}


int IncSatGC::external_get_fractional() {
    std::vector<Graph::VertexType> CliSAT_clique;
    namespace bp = boost::process::v1;

    //first check that external CliSAT binary for cliques exists
    std::filesystem::path EC_path;
#ifdef EC_BINARY_PATH
    EC_path = EC_BINARY_PATH;
    if (not std::filesystem::exists(EC_path)) {
        throw std::runtime_error("The external exactcolors binary does not exist at the specified path.");
    }
#else
    throw std::runtime_error("Path to exactcolors was not provided");
#endif

    std::string tmp_name;
    if(has_removed_vertices_in_reduction){
        //write current problem graph to file, so we can run EC
        tmp_name = "tmp_reduced_" + std::string(options.filename);
        graph.write_dimacs(tmp_name);
    }

    std::vector<std::string> args(2);
    args[0] = (has_removed_vertices_in_reduction ? tmp_name : options.filepath);
    args[1] = "-s 0";

    std::string cmd_output;
    bp::ipstream out_stream; // stream for command output

    auto start = Statistics::childCpuTime();
    // execute the command and redirect the output to out_stream
    bp::child c(EC_path.string(), args, bp::std_out > out_stream);

    auto TIMEOUT = std::chrono::seconds(initial_fractional_timeout);
    bool finished = c.wait_for(TIMEOUT);

    //add child process time to stats
    stats.add_fractional_time(Statistics::childCpuTime() - start);

    if (not finished) {
        c.terminate();
        if(has_removed_vertices_in_reduction){
            std::filesystem::remove(tmp_name);
        }
        flag_fractional_timed_out = true;
        return 0; //ec did not produce a lower bound
    }
    else {
        assert(stats.duration_of(Statistics::PreprocessingFractional).count() <= TIMEOUT.count());
    }

    std::regex  re{R"(LB\s+(\d+))", std::regex::icase};
    std::smatch m;
    int opt_colors = -1;

    std::string line;
    while (std::getline(out_stream, line)) {
        if (std::regex_search(line, m, re)) {
            opt_colors = std::stoi(m[1]);
        }
    }
    if (opt_colors < 0) {
        throw std::runtime_error("CliSAT finished but did not print 'Opt Colors:' line");
    }

    if(has_removed_vertices_in_reduction){
        //delete written temporary graph
        std::filesystem::remove(tmp_name);
    }
    return opt_colors;
}

void IncSatGC::notify_new_bound(const bool res, const int num_colors) {
    if(options.verbosity >= Options::Normal){
        std::cout << "Result: " << (res ? "Satisfiable" : "Unsatisfiable") << " for " << num_colors << " colors at "
                  << stats.current_total_time() << " and " << Statistics::memUsage() << " mb memory.\n";
    }
    //collect info of the algorithm up till this bound was found
    collect_bound_information(num_colors, res);
}

void IncSatGC::notify_already_used_fewer_colors(const int num_colors) {
    if(options.verbosity >= Options::Normal) {
        std::cout << "Result: Satisfiable for " << num_colors << " colors at "
                  << stats.current_total_time() << " (Model already used fewer colors).\n";
    }
    //collect info of the algorithm up till this upper bound was found
    collect_bound_information(num_colors, true);
}

void IncSatGC::notify_lower_bound(const int num_colors) {
    if(lower_bound < num_colors){
        lower_bound = num_colors;
        stats.lower_bound = num_colors;
        stats.tt_lower_bound = stats.current_total_time();
    }
}

void IncSatGC::notify_clique_lb(const int num_colors) {
    if(clique_lower_bound < num_colors){
        clique_lower_bound = num_colors;
        stats.clique_lower_bound = num_colors;
        notify_lower_bound(num_colors);
    }
}

void IncSatGC::notify_mycielsky_lb(const int num_colors) {
    if(mc_lower_bound < num_colors){
        mc_lower_bound = num_colors;
        stats.mc_lower_bound = num_colors;
        notify_lower_bound(num_colors);
    }
}

void IncSatGC::notify_fractional_lb(const int num_colors) {
    if(frac_lower_bound < num_colors){
        frac_lower_bound = num_colors;
        stats.frac_lower_bound = num_colors;
        notify_lower_bound(num_colors);
    }
}

void IncSatGC::notify_upper_bound(const int num_colors) {
    if(num_colors < upper_bound) {
        upper_bound = num_colors;
        stats.upper_bound = num_colors;
        stats.tt_upper_bound = stats.current_total_time();
    }
}

void IncSatGC::notify_heuristic_ub(const int num_colors) {
    if(num_colors < heuristic_bound) {
        heuristic_bound = num_colors;
        stats.heuristic_bound = num_colors;
        notify_upper_bound(num_colors);
    }
}

void IncSatGC::collect_bound_information(const int bound, const bool res) {
    //some bound information tracking
    stats.bound_information.push_back(stats.produce_bound_info(bound,res));
}


void IncSatGC::print_single_k_solved_by_upper_bound(const int num_colors) const {
    if(options.verbosity >= Options::Normal){
        std::cout << "Result: " << "Satisfiable for " << num_colors
                  << " colors since there is an upper bound of " << upper_bound << ".\n";
    }
}

void IncSatGC::print_single_k_solved_by_lower_bound(const int num_colors) const {
    if(options.verbosity >= Options::Normal){
        std::cout << "Result: " << "Unsatisfiable for " << num_colors
                  << " colors since there is a lower bound of " << lower_bound << ".\n";
    }
}

void IncSatGC::print_short_cegar_stats() const {
    if(options.verbosity >= Options::Verbose) {
        std::cout << "c Stats: found " << stats.num_total_conflicts << " conflicts in " << stats.num_cegar_iterations
                  << " CEGAR iterations, taking overall " << stats.duration_of(Statistics::CEGAR)
                  << " while SAT solver took " << stats.duration_of(Statistics::SatSolver) << "\n";
    }
}

void IncSatGC::print_sat_size() {
    if(options.verbosity >= Options::Debug){
        std::cout << "c Debug: " << get_num_vars() << " vars and " << get_num_clauses() << " clauses, "
        << Statistics::memUsage() << " mb memory.\n";
    }
    //mostly call print_sat_size when something changes so update those stats here
    //the number of variables and clauses should only increase, so this update should be sufficient
    stats.num_vars = get_num_vars();
    stats.num_clauses = get_num_clauses();
}


bool IncSatGC::is_valid_coloring(const ColorMap &colors) const {
    //iterate over all edges and make sure that their end vertices have different colors
    std::vector<Graph::VertexType> edge_list = graph.elist();
    for(int index =  0; 2*index < static_cast<int>(edge_list.size()); index++) {
        Graph::VertexType i = edge_list[2 * index];
        Graph::VertexType j = edge_list[2 * index + 1];
        if(colors[i] == colors[j]){
            std::cout << "Vertex " << i << " and " << j << " have the same color!\n";
            return false; //neighbors have same color
        }
    }
    //check that no vertex is uncolored
    for (int i = 0; i < num_vertices; ++i) {
        if(colors[i] == NoColor){
            std::cout << "Vertex " << i << " is uncolored!\n";
            return false;
        }
    }
    //is an optimal coloring
    int num_colors_used = *std::max_element(colors.begin(), colors.end()) + 1;
    if (num_colors_used > upper_bound){
        std::cout << "The coloring uses " << num_colors_used << " and upper bound is " << upper_bound << "\n";
        return false;
    }
    return true;
}


bool IncSatGC::is_valid_coloring(const Graph::Coloring &coloring) const {
    //every vertex is colored
    for(int v = 0; v < graph.ncount(); ++v){
        bool has_color = false;
        for (const Graph::ColorClass &cc : coloring) {
            if(cc.count(v)){
                has_color = true;
                break;
            }
        }
        if(not has_color){
            std::cout << "Vertex " << v << " is not assigned a color\n";
            return false;
        }
    }
    //every colorclass is an independent set
    auto neighborlist = graph.get_neighbor_list();
    for (const Graph::ColorClass &cc : coloring) {
        for (auto i = cc.begin(); i != cc.end(); ++i) {
            for (auto j = std::next(i); j != cc.end(); ++j) {
                //check that vertices in class are not adjacent
                if(std::find(neighborlist[*i].begin(), neighborlist[*i].end(), *j) != neighborlist[*i].end()){
                    std::cout << "Vertices " << *i << " and " << *j << " have the same color but are adjacent\n";
                    return false;
                }
            }
        }
    }
    //is an optimal coloring
    if (static_cast<int>(coloring.size()) != upper_bound){
        std::cout << "The coloring uses " << coloring.size() << " and upper bound is " << upper_bound << "\n";
        return false;
    }
    return true;
}


ColorMap IncSatGC::obtain_coloring_from_direct_encoding(int num_colors) const {
    //same method for direct and partial order encoding, the assignment variables are the same
    ColorMap colors(num_vertices, NoColor);
    Model model = get_model();
    assert(static_cast<int>(model.size()) == num_colors * num_vertices);
    if(options.encoding == Options::AssignmentEncoding) {
        //have k|V| variables x_v,j  for vertices 0<=v<=num_vertices-1 and color j
        //map x_v,j to variable with index k*v + j, first index is 0
        for (int v = 0; v < num_vertices; v++) {
            for (int j = 0; j < num_colors; j++) {
                int index = num_colors * v + j;
                if(model[index]){ //vertex i gets color j
                    colors[v] = j; //possibly overwrite color but that doesn't matter, simply choose one of the assigned
                }
            }
        }
    }
    else {
        assert(options.encoding == Options::PartialOrderEncoding);
        //extract coloring from solution to partial order encoding
        //vertex v has color i if y_v,i-1 is true and y_v,i is false
        for (int v = 0; v < num_vertices; v++) {
            if(not model[num_colors * v + 0]) { //special consideration for color 0
                colors[v] = 0;
                continue;
            }
            for (int j = 1; j < num_colors; j++) {
                int index = num_colors * v + j;
                if(model[index - 1] and not model[index]){
                    colors[v] = j;
                    break;
                }
            }
        }
    }
    assert(is_valid_coloring(colors));
    return colors;
}

int IncSatGC::get_num_colors_from_model() const {
    //if we obtain a model of the zykov encoding clauses, count the number of colors needed in the model
    // this might be smaller than the current at most k constraint, and we can skip a few iterations

    //look at all vertices and see if it has the same color as some previously used vertex
    //if yes, reduce number of needed colors, otherwise add it and all same colored vertices to the list
    std::unordered_set<int> same_colors;
    int num_colors = num_vertices;
    Model model = get_model();
    for(int i = 0; i < num_vertices; ++i){
        if(same_colors.count(i)) {//vertex already colored, do not need new color
            num_colors--;
        }
        for (int j : complement_graph_adjacency[i] ) {
            if(i > j){
                //only want to look at (i,j) once
                continue;
            }
            if(get_model_ij(model, i, j)){
                same_colors.insert(j);
            }
        }
    }
    return num_colors;
}


ColorMap IncSatGC::obtain_coloring_from_model() const {
    //look at all vertices and see if it has the same color as some previously used vertex
    //if yes, reduce number of needed colors, otherwise add it and all same colored vertices to the list
    int color_count = -1;
    ColorMap colors(num_vertices, NoColor);
    Model model = get_model();
    for(int i = 0; i < num_vertices; ++i){
        if(colors[i] == NoColor) {//vertex not colored, need a new color
            color_count++;
            colors[i] = color_count;//remember new color of vertex
        }
        for (int j : complement_graph_adjacency[i] ) {
            if(i > j){
                //only want to look at (i,j) once
                continue;
            }
            if(get_model_ij(model, i, j)){
                colors[j] = colors[i];//remember same color as i
            }
        }
    }
    assert(is_valid_coloring(colors));
    return colors;
}



void IncSatGC::write_optimal_coloring() {
    assert(not options.coloringfilepath.empty());
    //write found coloring to file
    std::string out_name = options.coloringfilepath;
    std::ofstream out_file = std::ofstream(out_name);
    if(!out_file.is_open()){
        std::cout << "Error: could not open file to write coloring to\n";
        return;
    }
    //coloring file header
    out_file << "c Coloring of graph " << options.filename << "\n"
             << "c Syntax: \n"
             << "c      p <problem name> <number of colors>\n"
             << "c      s <list of integers (= a color class)>\n";
    out_file << "p " << options.filename << " " << upper_bound << "\n";

    //convert colormap to collection of color classes
    Graph::Coloring coloring(upper_bound);
    if(not solved_in_preprocessing) {
        for (int i = 0; i < num_vertices; ++i) {
            //vertex is stored in its colorclass
            coloring[current_best_coloring[i]].insert(i);
        }
    }
    else {
        //graph was completely solved in reductions, or we computed (and proved) an optimal heuristic coloring
        coloring = heuristic_coloring;
    }
    //check we have a valid coloring
    assert(not coloring.empty());
    assert(std::all_of(coloring.begin(), coloring.end(),
        [this](const auto &cc){return not cc.empty();})
        or solved_in_preprocessing //can have empty colorclasses if solved entirely by reductions
        or options.strategy == Options::SingleK); //can have empty colorclasses if solving for a single (non-smallest) k

    //if we permuted the graph, reverse the permutation to find a coloring of the original graph
    if(options.use_clique_in_ordering and (not solved_in_preprocessing)){
        graph = graph.perm_graph(graph_vertex_ordering);
        for (int i = 0; i < static_cast<int>(coloring.size()); ++i) {
            //create copy to iterate over
            Graph::ColorClass cclass = coloring[i];
            //erase old element v, insert perm[v]
            // split up because we might insert a new vertex v before erasing old vertex v
            for (Graph::VertexType v : cclass) {
                coloring[i].erase(v);
            }
            for (Graph::VertexType v : cclass) {
                coloring[i].insert(graph_vertex_ordering[v]);
            }
            assert(cclass.size() == coloring[i].size());
        }
    }

    //vertices removed in preprocessing, insert them again and update + extend coloring
    if(has_removed_vertices_in_reduction){
        std::vector<Graph::VertexType> remapping = graph.get_remapping();
        std::vector<Graph::VertexType> inverse_remapping(remapping.size());
        for (int i = 0; i < static_cast<int>(remapping.size()); ++i) {
            if(remapping[i] != Graph::UndefVertex) {
                inverse_remapping[remapping[i]] = i;
            }
            else{
                inverse_remapping[i] = Graph::UndefVertex;
            }
        }
        //update coloring with inverse remapping
        //i.e. if original vertex v is sent to u with color c, erase u and assign v to color c
        for (int i = 0; i < static_cast<int>(coloring.size()); ++i) {
            //create copy to iterate over
            Graph::ColorClass cclass = coloring[i];
            //erase old element v, insert remapping[v]
            // split up because we might insert a new vertex v before erasing old vertex v
            for (Graph::VertexType v : cclass) {
                coloring[i].erase(v);
            }
            for (Graph::VertexType v : cclass) {
                assert(inverse_remapping[v] != Graph::UndefVertex);
                coloring[i].insert(inverse_remapping[v]);
            }
            assert(cclass.size() == coloring[i].size());
        }

        //now undo the reductions and get the set of recovered vertices in order of last removed to first
        std::vector<Graph::VertexType> recovered_vertices = graph.recover_reductions();

        //still need to find a color for the vertices recovered in reduction,
        // possible to color them without an additional color though
        auto neighborlist = graph.get_neighbor_list();
        for (Graph::VertexType uncolored : recovered_vertices) {
            bool assigned_color = false;
            for (Graph::ColorClass &cc : coloring) {
                if(not cc.count(uncolored)){
                    //also check that no adjacent vertex is in colorclass
                    bool can_be_inserted = true;
                    for (Graph::VertexType nb : neighborlist[uncolored]) {
                        if(cc.count(nb)){
                            can_be_inserted = false;
                            break;
                        }
                    }
                    if(can_be_inserted){
                        cc.insert(uncolored);
                        assigned_color = true;
                        break;
                    }
                }
            }
            assert(assigned_color);
        }
    }

    //debug/assert that we have a valid coloring
    if(not is_valid_coloring(coloring)) {
        throw std::runtime_error("Optimal coloring recovered after reductions was not valid.");
    }


    //finally write color classes
    for (int i = 0; i < upper_bound; ++i) {
        out_file << "s";
        for (int v : coloring[i]) {
            out_file << " " << v;
        }
        out_file << "\n";
    }
    out_file.close();
}


void IncSatGC::build_direct_encoding(int num_colors) {
    if (options.encoding == Options::AssignmentEncoding){
        build_assignment_encoding(num_colors);
    }
    else if(options.encoding == Options::PartialOrderEncoding){
        build_partial_order_encoding(num_colors);
    }
    else {
        throw std::runtime_error("Invalid encoding option in build_direct_encoding");
    }
}

bool IncSatGC::direct_encoding_single_k() {
    assert(options.encoding == Options::AssignmentEncoding or options.encoding == Options::PartialOrderEncoding);
    if (not options.specific_num_colors.has_value()) {
        throw std::runtime_error("Chose single k as search strategy but no k was given.");
    }
    int num_colors = options.specific_num_colors.value();
    if(num_colors >= upper_bound){
        print_single_k_solved_by_upper_bound(num_colors);
        return true;
    }
    if(lower_bound > num_colors){
        print_single_k_solved_by_lower_bound(num_colors);
        return false;
    }
    //build encoding and solve
    build_direct_encoding(num_colors);
    print_sat_size();
    if(options.write_cnf_only){
        write_cnf();
        write_and_cleanup();
        exit(0);
    }
    bool res = run_solver();
    notify_new_bound(res, num_colors);
    if(res){
        notify_upper_bound(num_colors);
        current_best_coloring = obtain_coloring_from_direct_encoding(num_colors);
    }
    else{
        notify_lower_bound(num_colors + 1);
    }
    return res;
}

int IncSatGC::direct_encoding_top_down() {
    assert(options.encoding == Options::AssignmentEncoding or options.encoding == Options::PartialOrderEncoding);
    if(options.specific_num_colors.has_value()){
        //a number of colors was given, use that as starting point bound in top-down approach
        notify_upper_bound(options.specific_num_colors.value() + 1);
    }
    //build large problem with upper_bound - 1 many colors
    int num_colors = upper_bound - 1;
    int full_num_colors = num_colors;
    build_direct_encoding(full_num_colors);
    print_sat_size();
    //start top down solving
    while(lower_bound != upper_bound){
        bool res = run_solver();
        notify_new_bound(res, num_colors);
        if (res){
            notify_upper_bound(num_colors);
            //found a k-coloring, decrease number of colors and continue
            //update best found coloring
            current_best_coloring = obtain_coloring_from_direct_encoding(full_num_colors);
            //check if model actually uses fewer colors already
            int num_colors_used = *std::max_element(current_best_coloring.begin(), current_best_coloring.end()) + 1;
            if(num_colors_used < num_colors){
                notify_upper_bound(num_colors_used);
                //update model according to actual num_colors_used
                notify_already_used_fewer_colors(num_colors_used);
                //constrain problem to not use variables with colors larger than num_colors_used, add them as unit literals
                for (int color = num_colors; color >= num_colors_used; --color) {
                    direct_encoding_top_down_fix_variables(full_num_colors, color);
                }
                num_colors = num_colors_used;
            }
            else{ //is equal to previous in case of num_colors_used == num_colors
                //constrain problem to not use variables with last color, add them as unit literals
                direct_encoding_top_down_fix_variables(full_num_colors, num_colors);
            }
            num_colors--;
        }
        else {//not satisfiable. chromatic number is current upper bound, stop loop
            notify_lower_bound(upper_bound);
        }
    }
    //chromatic number is lower_bound=upper_bound
    stats.solved = true;
    return upper_bound;
}

void IncSatGC::direct_encoding_top_down_fix_variables(const int full_num_colors, const int colors) {
    //fix the assignment variables x_v,colors-1 to false,
    if(options.encoding == Options::AssignmentEncoding) {
        for (int v = 0; v < num_vertices; ++v) {
            solver->addClause(mkLit(full_num_colors * v + (colors - 1), true));
        }
    }
    //in the partial order encoding, have to fix y_v,colors-2
    //since we want that no vertex can have colors larger than colors-2
    if(options.encoding == Options::PartialOrderEncoding){
        for (int v = 0; v < num_vertices; ++v) {
            solver->addClause(mkLit(full_num_colors * v + (colors - 2), true));
        }
    }
}

int IncSatGC::direct_encoding_bottom_up() {
    assert(options.encoding == Options::AssignmentEncoding or options.encoding == Options::PartialOrderEncoding);
    if(options.specific_num_colors.has_value()){
        //a number of colors was given, use that as starting point bound in bottom-up approach
        notify_lower_bound(options.specific_num_colors.value());
    }
    int num_colors = lower_bound;
    //start bottom up solving
    while ( lower_bound != upper_bound ) {
        //build and solve direct encoding for each k, there is no incremental procedure between different colors
        reset_SAT_solver();
        build_direct_encoding(num_colors);
        print_sat_size();
        bool res = run_solver();
        notify_new_bound(res, num_colors);
        if (res){
            //found a k-coloring, became satisfiable
            notify_upper_bound(num_colors); //num_colors == lower_bound == first satisfiable k
            //update best found coloring (will be an optimal one)
            current_best_coloring = obtain_coloring_from_direct_encoding(num_colors);
            assert(*std::max_element(current_best_coloring.begin(), current_best_coloring.end()) + 1 == num_colors);
            assert(lower_bound == upper_bound);
        }
        else {//not satisfiable, increase number of colors and continue
            notify_lower_bound(num_colors + 1); //+1 because num_colors was UNSAT
            num_colors++;
        }
    }
    //chromatic number is lower_bound=upper_bound
    stats.solved = true;
    return upper_bound;
}

int IncSatGC::do_assignment_encoding() {
    switch (options.strategy) {
        case Options::SingleK:
            assignment_encoding_single_k();
            break;
        case Options::TopDown:
            assignment_encoding_top_down();
            break;
        case Options::BottomUp:
            assignment_encoding_bottom_up();
            break;
        default:
            throw std::runtime_error("Invalid option for search strategy.");
    }
    return 0;
}


void IncSatGC::build_assignment_encoding(const int num_colors) {
    stats.start_phase(Statistics::BuildEncoding);
    int num_vars = num_colors * num_vertices;
    //have k|V| variables x_v,j  for vertices 0<=v<=num_vertices-1 and color j
    //map x_v,j to variable with index k*v + j, first index is 0
    add_vars(num_vars);

    //found a clique in preprocessing, use it to fix some variables by assigning clique vertices different colors
    if(not options.disable_preprocessing) {
        for (int index = 0; index < static_cast<int>(clique.size()); ++index) {
            for (int color = 0; color < num_colors; ++color) {
                if(index == color){
                    //set x_clique[index],color to true
                    solver->addClause(mkLit( num_colors * (clique[index]) + color, false));
                }
                else{
                    //set x_clique[index],color to false
                    solver->addClause(mkLit( num_colors * (clique[index]) + color, true));
                }
            }
        }
    }

    assert(not tmp_clause.size());
    //add at-least-one-color constraints: for each v, OR_j x_v,c
    for (int v = 0; v < num_vertices; v++) {
        for (int c = 0; c < num_colors; c++) {
            tmp_clause.push(mkLit(num_colors * v + c, false));
        }
        solver->addClause(tmp_clause);
        tmp_clause.clear();
    }

    //add different-colors constraint
    //for edge (u,v) and color c, ensure that one of x_u,c or x_v,c is false
    std::vector<Graph::VertexType> edge_list = graph.elist();
    for (int c = 0; c < num_colors; c++) {
        for (int index = 0; 2 * index < static_cast<int>(edge_list.size()); index++) {
            tmp_clause.push(mkLit(num_colors * (edge_list[2 * index]) + c, true));
            tmp_clause.push(mkLit(num_colors * (edge_list[2 * index + 1]) + c, true));
            solver->addClause(tmp_clause);
            tmp_clause.clear();
        }
    }

    //optionally add at-most-one color constraint (-x_v,i or -x_v,j)
    if(options.assignment_encoding_amo) {
        for (int v = 0; v < num_vertices; v++) {
            for (int i = 0; i < num_colors; i++) {
                for (int j = i + 1; j < num_colors; j++) {
                    tmp_clause.push(mkLit(num_colors * v + i, true));
                    tmp_clause.push(mkLit(num_colors * v + j, true));
                    solver->addClause(tmp_clause);
                    tmp_clause.clear();
                }
            }
        }
    }

    if (options.use_clique_in_ordering){ //only applicable if vertex ordering starts with clique
        //also add simple symmetry breaking clauses: -x_v,i for v < i
        for (int i = 1; i < num_colors; i++) {
            for (int v = 0; v < i; v++) {
                if (std::find(clique.begin(), clique.end(), v) != clique.end()) {
                    continue; //already fixed vertices in clique, ignore them here
                }
                solver->addClause(mkLit(num_colors * v + i, true));
            }
        }
    }


    //could also consider symmetry breaking clauses as in Van Gelder or Van Hoeve paper
    stats.end_phase(Statistics::BuildEncoding);
}


void IncSatGC::build_partial_order_encoding(const int num_colors) {
    stats.start_phase(Statistics::BuildEncoding);
    int num_vars = num_colors * num_vertices;
    //have k|V| variables y_v,i  for vertices 0<=v<=num_vertices-1 and color i
    //y_v,i means vertex v has color greater than i, -y_v,i means color at most i
    //map y_v,i to variable with index k * v + i, first index is 0
    add_vars(num_vars);
    tmp_clause.clear();
    std::vector<Graph::VertexType> edge_list = graph.elist();


    //fix assignments for vertices in clique
    //for ith vertex v in clique, set y_v,j true for j<i and false for j>= i
    for (int index = 0; index < static_cast<int>(clique.size()); ++index) {
        for (int j = 0; j < index; j++) {
            //set y_clique[index],j to true
            solver->addClause(mkLit( num_colors * (clique[index]) + j, false));
        }
        for (int j = index; j < num_colors; j++) {
            //set y_clique[index],j to false
            solver->addClause(mkLit( num_colors * (clique[index]) + j, true));
        }
    }

    //add unit clause -y,v,k for largest color k
    for (int v = 0; v < num_vertices; v++) {
        solver->addClause(mkLit( num_colors * v + (num_colors - 1), true));
    }

    //ensure that v has color smaller k means it also has color smaller k+1
    for (int v = 0; v < num_vertices; v++) {
        for (int c = 0; c < num_colors - 1; c++) {
            tmp_clause.push(mkLit(num_colors * v + c, false));
            tmp_clause.push(mkLit(num_colors * v + (c + 1), true));
            solver->addClause(tmp_clause);
            tmp_clause.clear();
        }
    }

    //ensure that adjacent vertices get different colors
    for (int index = 0; 2 * index < static_cast<int>(edge_list.size()); index++) {
        int u = edge_list[2 * index];
        int v = edge_list[2 * index + 1];
        //for color 0
        tmp_clause.push(mkLit(num_colors * u + 0, false));
        tmp_clause.push(mkLit(num_colors * v + 0, false));
        solver->addClause(tmp_clause);
        tmp_clause.clear();
        //for other colors
        for (int c = 0; c < num_colors - 1; c++) {
            tmp_clause.push(mkLit(num_colors * u + c, true));       //-y_u,c
            tmp_clause.push(mkLit(num_colors * u + (c + 1), false));//y_u,c+1
            tmp_clause.push(mkLit(num_colors * v + c, true));       //-y_v,c
            tmp_clause.push(mkLit(num_colors * v + (c + 1), false));//y_v,c+1
            solver->addClause(tmp_clause);
            tmp_clause.clear();
        }
    }

    if (options.use_clique_in_ordering) {//only applicable if vertex ordering starts with clique
        //also add simple symmetry breaking clauses: -y_i,i for i = 1...k
        for (int i = 0; i < num_colors; i++) {
            if(std::find(clique.begin(), clique.end(), i) != clique.end()) {
                continue; //already fixed vertices in clique, ignore them here
            }
            solver->addClause(mkLit( num_colors * i + i, true));
        }
    }
    stats.end_phase(Statistics::BuildEncoding);
}


bool IncSatGC::assignment_encoding_single_k() {
    return direct_encoding_single_k();
}

int IncSatGC::assignment_encoding_top_down(){
    return direct_encoding_top_down();
}

int IncSatGC::assignment_encoding_bottom_up() {
    return direct_encoding_bottom_up();
}

int IncSatGC::do_partial_order_encoding() {
    switch (options.strategy) {
        case Options::SingleK:
            partial_order_encoding_single_k();
            break;
        case Options::TopDown:
            partial_order_encoding_top_down();
            break;
        case Options::BottomUp:
            partial_order_encoding_bottom_up();
            break;
        default:
            throw std::runtime_error("Invalid option for search strategy.");
    }
    return 0;
}

bool IncSatGC::partial_order_encoding_single_k() {
    return direct_encoding_single_k();
}

int IncSatGC::partial_order_encoding_top_down() {
    return direct_encoding_top_down();
}

int IncSatGC::partial_order_encoding_bottom_up() {
    return direct_encoding_bottom_up();
}


void IncSatGC::write_full_maxsat_encoding() {
    stats.start_phase(Statistics::BuildEncoding);
    //write wcnf to file
    std::string out_name = options.filename + "_MaxSAT.wcnf";
    wcnf = std::ofstream(out_name);
    if(!wcnf.is_open()){
        throw std::runtime_error("Could not open file to write wcnf.");
    }
    //initialise indices of s_ij and c_j in member fields
    initialise_variable_indices(1);//wncf format starts with 1

    //initialising and indices done, start writing formula
    wcnf << "c wcnf formula that uses the zykov encoding and all transitivity constraints \n"
         << "c to compute the chromatic number with this MaxSAT formula.\n";
    //compute header, "p wcnf nbvars nbclauses top" where top marks hard clauses
    int nbvars = c_indices.back();
    int nbclauses = -1; //can't say yet how many clauses, seems to be fine for solvers though
    int top = graph.ncount();
    wcnf << "p wcnf " << nbvars << " " << nbclauses << " " << top << "\n";

    int count_clauses = 0;// num_clauses is counting too many clauses because we get rid of a lot already

    //don't need to add unit literals as we do not add variables for edges, only for non-edges

    //now add all transitivity constraints, that is for i<j<k add
    // transitivity(i,j,k) & transitivity(j,i,k) & transitivity(i,k,j)
    write_all_transitivity(count_clauses);
    //next are the definitions of c_j as in equation 2
    write_cj_definition(count_clauses);

    //for MaxSAT, add the c_j as the soft clauses
    for (int j = 1; j < num_vertices; ++j) {
        wcnf << 1 << " " << -c_indices[j - 1] << " 0\n";
        count_clauses++;
    }
    // no need for at most k-1 constraint since we solve the MaxSAT problem

    if(options.verbosity >= Options::Debug) {
        std::cout << "c Debug: " << nbvars << " vars and " << count_clauses << " clauses\n";
    }
    std::cout << "c Debug: Wrote MaxSAT formula to " << out_name << "\n";
    wcnf.close();
    stats.end_phase(Statistics::BuildEncoding);
}

void IncSatGC::write_all_transitivity(int &count_clauses) {
    for (int i = 0; i < num_vertices; ++i) {
        for (int j = i+1; j < num_vertices; ++j) {
            for (int k = j+1; k < num_vertices; ++k) {
                //first case transitivity(i,j,k)
                write_transitivity(count_clauses, i, j, k);
                //second case transitivity(j,i,k)
                write_transitivity(count_clauses, j, i, k);
                //third case transitivity(i,k,j)
                write_transitivity(count_clauses, i, k, j);
            }
        }
    }
}


void IncSatGC::write_transitivity(int &count_clauses, int i, int j, int k) {
    int top = sij_indices.dimension; // top is num vertices is dimension
    if((sij_indices.get(i, j) != -1) and (sij_indices.get(j, k)) != -1){
        //both i,j and j,k are not an edge, need to ensure transitivity i,j,k
        //unless i,k is an edge, then the clause changes slightly
        wcnf << top << " " << -sij_indices.get(i, j) << " " << -sij_indices.get(j, k) << " ";
        if(sij_indices.get(i, k) != -1){ //i,k also a non-edge, include in clause
            wcnf << sij_indices.get(i, k) << " ";
        }
        wcnf << "0 \n";
        count_clauses++;
    }
}

void IncSatGC::write_cj_definition(int &count_clauses) {
    int top = num_vertices;
    for (int j = 1; j < num_vertices ; ++j) {
        //for each j, add clause c_j <-> AND -s_ij
        //which in cnf is AND(-c_j v -s_ij) and (ORs_ij v c_j)
        //do first clauses while collecting literals for second
        std::vector<int> second_clause;
        second_clause.reserve(j);
        for (int i = 0; i < j; ++i) {
            if(sij_indices.get(i, j) == -1){
                //is an edge, thus leave it out
                continue;
            }
            wcnf << top << " " << -c_indices[j - 1] << " " << -sij_indices.get(i, j) << " 0 \n";
            count_clauses++;

            second_clause.push_back(sij_indices.get(i, j));
        }
        second_clause.push_back(c_indices[j - 1]);
        wcnf << top << " ";
        for (int lit : second_clause) {
            wcnf << lit << " ";
        }
        wcnf << "0\n";
        count_clauses++;
    }
}


void IncSatGC::initialise_variable_indices(int start_index) {
    //we have variables s_i,j for all i,j in V, only consider i<j so that is n(n-1)/2 in upper triangular matrix
    //only look at i,j such that (i,j) is not an edge
    //and variables c_j for 1 <= j <= n-1. keep track of indices in vectors
    sij_indices = UpperTriangle(num_vertices);
    c_indices.resize(num_vertices - 1);
    int var_index = start_index;

    complement_graph_adjacency = graph.get_complement().get_neighbor_list();
    for (int i = 0; i < num_vertices; ++i) {
        for (int j : complement_graph_adjacency[i] ) {
            if(i > j){
                //only want to look at (i,j) once
                continue;
            }
            //i and j are not adjacent in original graph
            sij_indices.set(i,j, var_index);
            var_index++;
        }
    }
    stats.num_sij_vars = var_index - start_index;

    for (int j = 1; j < num_vertices; ++j) {
        c_indices[j - 1] = var_index;
        var_index++;
    }
    stats.num_cj_vars = num_vertices - 1;
    assert(stats.num_cj_vars == var_index - start_index - stats.num_sij_vars);
}

void IncSatGC::add_transitivity(int i, int j, int k) {
    assert(not tmp_clause.size());
    if((sij_indices.get(i, j) != -1) and (sij_indices.get(j, k)) != -1){
        //both i,j and j,k are not an edge, need to ensure transitivity i,j,k
        //unless i,k is an edge, then the clause changes slightly
        tmp_clause.push(mkLit(sij_indices.get(i, j), true));
        tmp_clause.push(mkLit(sij_indices.get(j, k), true));
        if(sij_indices.get(i, k) != -1){ //i,k also a non-edge, include in clause
            tmp_clause.push(mkLit(sij_indices.get(i, k), false));
        }
        solver->addClause(tmp_clause);
        tmp_clause.clear();
        stats.num_transitivity_clauses++;
    }
}



void IncSatGC::add_all_transitivity() {
    for (int i = 0; i < num_vertices; ++i) {
        for (int j = i+1; j < num_vertices; ++j) {
            for (int k = j+1; k < num_vertices; ++k) {
                //first case transitivity(i,j,k)
                add_transitivity(i, j, k);
                //second case transitivity(j,i,k)
                add_transitivity(j, i, k);
                //third case transitivity(i,k,j)
                add_transitivity(i, k, j);
            }
        }
    }
}

void IncSatGC::add_cj_definition() {
    //catch case that cardinality constraints are handled by cliques in zykov propagator
    if(options.encoding == Options::ZykovPropagator and options.disable_cardinality_constraints) {
        assert(options.use_clique_explanation_clauses);
        if(options.verbosity >= Options::Debug){
            std::cout << "c Debug: No cj variables were added\n";
        }
        return;
    }
    assert(not tmp_clause.size());
    vec<Lit> second_clause;
    num_removable_cj = 0;
    is_removable_cj = Bitset(num_vertices - 1);
    for (int j = 1; j < num_vertices ; ++j) {
        //for each j, add clause c_j <-> AND -s_ij
        //which in cnf is AND(-c_j v -s_ij) and (ORs_ij v c_j)
        //do first clauses while collecting literals for second
        second_clause.capacity(j);
        for (int i = 0; i < j; ++i) {
            if(sij_indices.get(i, j) == -1){
                //is an edge, thus leave it out
                continue;
            }
            tmp_clause.push(mkLit(sij_indices.get(i, j), true));
            tmp_clause.push(mkLit(c_indices[j - 1], true));
            solver->addClause(tmp_clause);
            tmp_clause.clear();

            second_clause.push(mkLit(sij_indices.get(i, j), false));
        }
        if(second_clause.size() == 0 and options.remove_trivial_cj) {
            //remember that c_j can be removed in cardinality constraints later as it must be set to true
            num_removable_cj++;
            is_removable_cj[j - 1] = true;
        }
        second_clause.push(mkLit(c_indices[j - 1], false));
        solver->addClause(second_clause);
        second_clause.clear();

    }
    if(options.verbosity >= Options::Debug and options.remove_trivial_cj){
        std::cout << "c Debug: Found " << num_removable_cj << " c_j to be removable. [" << is_removable_cj << "]\n";
    }
    notify_lower_bound(num_removable_cj + 1);
}

void IncSatGC::add_zykov_encoding() {
    stats.start_phase(Statistics::BuildEncoding);
    //initialise indices of s_ij and c_j in member fields
    initialise_variable_indices(0);
    add_vars(stats.num_sij_vars + stats.num_cj_vars);

    //don't need to add unit literals as we do not add variables for edges, only for non-edges

    if (options.encoding == Options::FullEncoding) {
        //now add all transitivity constraints, that is for i<j<k add
        //transitivity(i,j,k) & transitivity(j,i,k) & transitivity(i,k,j)
        add_all_transitivity();
    }

    //next are the definitions of c_j as in equation 2
    add_cj_definition();

    //if we are using the external propagator for solving, take care of initialising it
    if(options.encoding == Options::ZykovPropagator) {
        init_zykov_propagator();
    }
    stats.end_phase(Statistics::BuildEncoding);
}

void IncSatGC::add_at_most_k(int k) {
    assert(options.strategy == Options::SingleK or options.non_incremental_solving);
    assert(not cj_literals.size());
    //catch case that cardinality constraints are handled by cliques in zykov propagator
    if(options.encoding == Options::ZykovPropagator and options.disable_cardinality_constraints) {
        assert(zykov_propagator != nullptr);
        zykov_propagator->update_num_colors(k + 1);
        assert(options.use_clique_explanation_clauses);
        if(options.verbosity >= Options::Debug){
            std::cout << "c Debug: No cardinality constraints added, coloring size is asserted by cliques\n";
        }
        return;
    }

    stats.start_phase(Statistics::BuildAtMostK);
    int num_vars = get_num_vars();
    int num_clauses = get_num_clauses();

    for (int j = 1; j < num_vertices; ++j) {
        if(options.remove_trivial_cj and is_removable_cj[j - 1]){
            continue; //do not add c_j to cardinality constraints
        }
        cj_literals.push(NSPACE::mkLit(c_indices[j - 1], false));
    }

    int at_most_k = k;
    if(options.remove_trivial_cj){
        at_most_k -= num_removable_cj;
        if(options.verbosity >= Options::Debug){
            std::cout << "c Debug: Removed " << num_removable_cj << " c_j from cardinality constraints\n";
        }
        stats.num_removable_cj = num_removable_cj;
        assert(static_cast<int>(cj_literals.size()) == (num_vertices - 1 - num_removable_cj));
    }

    encoder.encodeCardinality(solver.get(), cj_literals, at_most_k);
    if (not encoder.hasCardEncoding()){
        if(options.verbosity >= Options::Debug){
            if(at_most_k == 0){
                std::cout << "c Encoding: removed enough vertices that rhs became 0, added -c_j as unit clause for the rest\n";
            }
            if(at_most_k == cj_literals.size()){
                std::cout << "c Encoding: have an at most k out of k constraint, nothing to be done\n";
            }
        }
        assert(at_most_k == 0 or at_most_k == static_cast<int>(cj_literals.size()));
    }
    if(options.verbosity >= Options::Debug){
        std::cout << "c Encoding: additional " << (get_num_vars() - num_vars) << " vars and " << (get_num_clauses() - num_clauses)
                  << " clauses for encoding of at most " << at_most_k << " out of " << (cj_literals.size()) << "\n";
    }
    //only encode cardinality once, can clear the literal vector
    cj_literals.clear();

    //for propagator, set num_colors so it has that number too
    if(options.encoding == Options::ZykovPropagator) {
        assert(zykov_propagator != nullptr);
        zykov_propagator->update_num_colors(k + 1);
    }

    stats.num_vars_at_most_k += (get_num_vars() - num_vars);
    stats.num_clauses_at_most_k += (get_num_clauses() - num_clauses);
    stats.end_phase(Statistics::BuildAtMostK);
}


void IncSatGC::add_incremental_at_most_k(int k) {
    assert(options.strategy != Options::SingleK);
    //catch case that cardinality constraints are handled by cliques in zykov propagator
    if(options.encoding == Options::ZykovPropagator and options.disable_cardinality_constraints) {
        assert(zykov_propagator != nullptr);
        zykov_propagator->update_num_colors(k + 1);
        assert(options.use_clique_explanation_clauses);
        if(options.verbosity >= Options::Debug){
            std::cout << "c Debug: No cardinality constraints added, coloring size is asserted by cliques\n";
        }
        return;
    }

    stats.start_phase(Statistics::BuildAtMostK);
    //weird behaviours of cardinality encoding, if rhs = 0 unit clauses are added which is a problem for the BottomUp approach
    //this ensures that in BottomUp rhs is at least 1 and actual encoding is built instead of unit clauses
    bool need_extra_consideration = (options.strategy == Options::BottomUp and k == num_removable_cj);
    if(options.remove_trivial_cj and need_extra_consideration and cj_literals.size() == 0) {
        num_removable_cj--;
        if(options.verbosity >= Options::Debug){
            std::cout << "c Debug: had to do extra considerations for BottomUp, kept 1 cj literal that was removable\n";
        }
    }
    int num_vars = get_num_vars();
    int num_clauses = get_num_clauses();
    // careful, this breaks the normal encodeCardinality
    encoder.setIncremental(openwbo::_INCREMENTAL_ITERATIVE_);
    int at_most_k = k;
    if(cj_literals.size() == 0){ //have not yet created list of literals for constraint
        for (int j = 1; j < num_vertices; ++j) {
            if(options.remove_trivial_cj and is_removable_cj[j - 1] and not need_extra_consideration){
                continue; //do not add c_j to cardinality constraints
            }
            need_extra_consideration = false;
            cj_literals.push(NSPACE::mkLit(c_indices[j - 1], false));
        }
    }
    if(options.remove_trivial_cj){
        at_most_k -= num_removable_cj;
        stats.num_removable_cj = num_removable_cj;
        assert(static_cast<int>(cj_literals.size()) == (num_vertices - 1 - num_removable_cj));
    }

    if(not encoder.hasCardEncoding()){//build an encoding if not present
        if(options.verbosity >= Options::Debug){
            std::cout << "c Debug: no encoding present, try to create one\n";
        }
        assert(at_most_k >= 0 and (at_most_k <= static_cast<int>(cj_literals.size())));
        encoder.buildCardinality(solver.get(), cj_literals, at_most_k);
    }
    //can happen that we do not build an encoding because the constraint is trivial
    if(encoder.hasCardEncoding()) {//encoding present, update rhs
        encoder.incUpdateCardinality(solver.get(), cj_literals, at_most_k, encoder_assumptions);
    }

    if(options.verbosity >= Options::Debug){
        std::cout << "c Encoding: additional " << (get_num_vars() - num_vars) << " vars and " << (get_num_clauses() - num_clauses)
                  << " clauses for encoding of at most " << at_most_k << " out of " << (cj_literals.size()) << "\n";
    }

    //for propagator, set num_colors so it has that number too
    if(options.encoding == Options::ZykovPropagator) {
        assert(zykov_propagator != nullptr);
        zykov_propagator->update_num_colors(k + 1);
    }

    stats.num_vars_at_most_k += (get_num_vars() - num_vars);
    stats.num_clauses_at_most_k += (get_num_clauses() - num_clauses);
    stats.end_phase(Statistics::BuildAtMostK);
}

bool IncSatGC::zykov_encoding_run_solver() {
    if( options.encoding == Options::CEGAR) {
        //additional inner loop, cegar approach of adding new clauses until no more conflicts
        return search_add_conflicts_and_solve();
    }
    assert(options.encoding == Options::FullEncoding or options.encoding == Options::ZykovPropagator);
    //simply solve full model or let propagator and callbacks run
    return run_solver();
}


bool IncSatGC::zykov_encoding_single_k() {
    assert(options.encoding == Options::FullEncoding or options.encoding == Options::CEGAR or options.encoding == Options::ZykovPropagator);
    if (not options.specific_num_colors.has_value()) {
        throw std::runtime_error("Chose single k as search strategy but no k was given.");
    }
    int num_colors = options.specific_num_colors.value();
    if(num_colors >= upper_bound){
        print_single_k_solved_by_upper_bound(num_colors);
        return true;
    }
    if(lower_bound > num_colors){
        print_single_k_solved_by_lower_bound(num_colors);
        return false;
    }
    //initialise variable indices and add constraints and clauses defining the c_j
    // either add all transitivity constraints or NOT, done via the settings
    // or, initalise the external propagator
    add_zykov_encoding();
    //add the at most k-1 constraints on the c_j
    add_at_most_k(num_colors - 1);
    print_sat_size();
    if(options.write_cnf_only){
        write_cnf();
        write_and_cleanup();
        exit(0);
    }
    bool res = zykov_encoding_run_solver();
    notify_new_bound(res, num_colors);
    if(res){
        notify_upper_bound(num_colors);
        current_best_coloring = obtain_coloring_from_model();
    }
    else{
        notify_lower_bound(num_colors + 1);
    }
    return res;
}

int IncSatGC::zykov_encoding_top_down() {
    assert(options.encoding == Options::FullEncoding or options.encoding == Options::CEGAR or options.encoding == Options::ZykovPropagator);
    if(options.specific_num_colors.has_value()){
        //a number of colors was given, use that as starting point bound in top-down approach
        notify_upper_bound(options.specific_num_colors.value() + 1);
    }

    //initialise variable indices, add all the transitivity constraints and clauses defining the c_j, or zykov propagator
    add_zykov_encoding();
    //for at most k constraints, start with one less than the best known SAT result
    int num_colors = upper_bound - 1;
    //start top down solving
    while(lower_bound != upper_bound) {//outer loop, decreasing num colors
        add_incremental_at_most_k(num_colors - 1);
        print_sat_size();

        bool res = zykov_encoding_run_solver();
        notify_new_bound(res, num_colors);
        if (res){
            notify_upper_bound(num_colors);
            //found a k-coloring, decrease number of colors and continue
            //check if model actually uses fewer colors already
            current_best_coloring = obtain_coloring_from_model();
            int num_colors_used = *std::max_element(current_best_coloring.begin(), current_best_coloring.end()) + 1;
            if(num_colors_used < num_colors){
                notify_upper_bound(num_colors_used);
                notify_already_used_fewer_colors(num_colors_used);
                num_colors = num_colors_used;
            }
            num_colors--;
        }
        else {//not satisfiable. chromatic number is current upper bound, stop loop
            notify_lower_bound(upper_bound);
        }
    }
    //chromatic number is lower_bound=upper_bound
    stats.solved = true;
    return upper_bound;
}

int IncSatGC::zykov_encoding_bottom_up() {
    assert(options.encoding == Options::FullEncoding or options.encoding == Options::CEGAR or options.encoding == Options::ZykovPropagator);
    if(options.specific_num_colors.has_value()){
        //a number of colors was given, use that as starting point bound in bottom-up approach
        notify_lower_bound(options.specific_num_colors.value());
    }

    //initialise variable indices, add all the transitivity constraints and clauses defining the c_j, or zykov propagator
    add_zykov_encoding();
    int num_colors = lower_bound;
    //start bottom up solving
    while(lower_bound != upper_bound) {//outer loop, increasing num colors
        if (options.non_incremental_solving == true) {
            //disable incremental solving so rebuild encoding with each call
            reset_SAT_solver();
            add_zykov_encoding();
            add_at_most_k(num_colors - 1);
        }
        else {//do incremental solving so update cardinality constraints incrementally
            add_incremental_at_most_k(num_colors - 1);
        }

        print_sat_size();

        bool res = zykov_encoding_run_solver();
        notify_new_bound(res, num_colors);
        if (res){//found a k-coloring, became satisfiable
            notify_upper_bound(num_colors);
            current_best_coloring = obtain_coloring_from_model();
            assert(*std::max_element(current_best_coloring.begin(), current_best_coloring.end()) + 1 == num_colors);
            assert(lower_bound == upper_bound);
            //loop stops here
        }
        else {//not satisfiable, increase number of colors and continue
            notify_lower_bound(num_colors + 1);
            num_colors++;
        }
    }
    //chromatic number is lower_bound=upper_bound
    stats.solved = true;
    return upper_bound;
}



int IncSatGC::do_full_encoding() {
    switch (options.strategy) {
        case Options::SingleK:
            full_encoding_single_k();
            break;
        case Options::TopDown:
            full_encoding_top_down();
            break;
        case Options::BottomUp:
            full_encoding_bottom_up();
            break;
        default:
            throw std::runtime_error("Invalid option for search strategy.");
    }
    return 0;
}


bool IncSatGC::full_encoding_single_k() {
    return zykov_encoding_single_k();
}

int IncSatGC::full_encoding_top_down(){
    return zykov_encoding_top_down();
}

int IncSatGC::full_encoding_bottom_up(){
    return zykov_encoding_bottom_up();
}



int IncSatGC::do_cegar_approach() {
    switch (options.strategy) {
        case Options::SingleK:
            cegar_approach_single_k();
            break;
        case Options::TopDown:
            cegar_approach_top_down();
            break;
        case Options::BottomUp:
            cegar_approach_bottom_up();
            break;
        default:
            throw std::runtime_error("Invalid option for search strategy.");
    }
    return 0;
}

bool IncSatGC::cegar_approach_single_k() {
    return zykov_encoding_single_k();
}

int IncSatGC::cegar_approach_top_down() {
    return zykov_encoding_top_down();
}

int IncSatGC::cegar_approach_bottom_up() {
    return zykov_encoding_bottom_up();
}


bool IncSatGC::search_add_conflicts_and_solve() {
    ConflictList conflicts;
    bool res;
    do {
        conflicts.clear();
        res = run_solver();
        if(not res){ //instance UNSAT, no need to look for further conflicts
            break;
        }

        //check for and add conflict clauses
        stats.start_phase(Statistics::CEGAR);
        conflicts = find_conflicts();
        add_conflict_clauses(conflicts);
        stats.num_cegar_iterations++;
        stats.num_total_conflicts += static_cast<int>(conflicts.size());
        stats.store_num_conflicts.push_back(static_cast<int>(conflicts.size()));
        stats.end_phase(Statistics::CEGAR);

        if(stats.num_cegar_iterations % 100 == 0){
            print_short_cegar_stats();
        }
    } while ( not conflicts.empty());

    print_short_cegar_stats();

    return res;
}


ConflictList IncSatGC::find_conflicts() const {
    switch (options.checker) {
        case Options::NaiveChecker:
            return find_conflicts_naive();
        case Options::SparseTrianglesChecker:
            return find_conflicts_sparse_triangles();
        case Options::AllTrianglesChecker:
            return find_conflicts_all_triangles();
        case Options::PaperChecker:
            return find_conflicts_paper();
        default:
            throw std::runtime_error("Invalid option set for checker algorithm.");
    }
}

ConflictList IncSatGC::find_conflicts_naive() const {
    //simply tries all triples and checks for violated constraints, adds those triples
    Model model = get_model();
    ConflictList conflicts;
    for (int i = 0; i < num_vertices; ++i) {
        for (int j = i+1; j < num_vertices; ++j) {
            for (int k = j+1; k < num_vertices; ++k) {
                bool model_ij = get_model_ij(model, i, j);
                bool model_jk = get_model_ij(model, j, k);
                bool model_ik = get_model_ij(model, i, k);
                //first case transitivity(i,j,k)
                if( model_ij and model_jk and !model_ik){
                    conflicts.emplace_back(i,j,k);
                }
                //second case transitivity(j,i,k)
                if( model_ij and model_ik and !model_jk){
                    conflicts.emplace_back(j,i,k);
                }
                //third case transitivity(i,k,j)
                if( model_ik and model_jk and !model_ij){
                    conflicts.emplace_back(i,k,j);
                }
            }
        }
    }
    return conflicts;
}

Graph::NeighborList IncSatGC::get_pairing_graph(const Model &model) const {
    Graph::NeighborList pairing_graph(num_vertices);
    //add edge for each s_ij set to true in model, iterate over each non-edge to check
    for (int i = 0; i < num_vertices; ++i) {
        for (int j : complement_graph_adjacency[i] ) {
            if(i > j){
                //only want to look at (i,j) once
                continue;
            }
            if(get_model_ij(model, i, j)) {
                pairing_graph[i].push_back(j);
                pairing_graph[j].push_back(i);
            }
        }
    }
    return pairing_graph;
}

ConflictList IncSatGC::find_conflicts_sparse_triangles() const {
    //idea for this checker: same as find_conflicts_all_triangles but only add one conflcit for each visited vertex
    ConflictList conflicts;
    Model model = get_model();
    Graph::NeighborList pairing_graph = get_pairing_graph(model);

    //constructed graph with edges all the s_i set to true,
    // iterate over all its vertices
    // and for each pair of edges, check whether there is also the third edge in the triangle -> if not, conflict
    // runtime n max_d^2, where max degree should be relatively small, since we look only at all s_ij that are true
    for (int vertex_j = 0; vertex_j < num_vertices; vertex_j++) {
        //next, iterate over pairs of neighbours
        bool found_conflict = false;
        for(auto i_it = pairing_graph[vertex_j].begin(); i_it != pairing_graph[vertex_j].end() and not found_conflict; ++i_it){
            for (auto k_it = std::next(i_it); k_it != pairing_graph[vertex_j].end() and not found_conflict; ++k_it) {
                //from the pairing graph we know that s_ij and s_jk are true, check for s_ik
                if(not get_model_ij(model, *i_it, *k_it)){
                    //found a conflict: s_ij and s_jk are true but s_ik is not
                    conflicts.emplace_back(*i_it, vertex_j, *k_it);
                    //found conflict for vertex_j, continue with next vertex by exiting both inner loops
                    found_conflict = true;
                }
            }
        }
    }
    if(conflicts.empty()){
        assert(find_conflicts_naive().empty());
    }
    return conflicts;
}

ConflictList IncSatGC::find_conflicts_all_triangles() const {
    //idea for this checker: build new graph with s_ij as edges
    //for all pairs of edges, add their transitivity constraints
    //this covers all current conflicts
    //Basically the same as Connected component checker,
    ConflictList conflicts;
    Model model = get_model();
    Graph::NeighborList pairing_graph = get_pairing_graph(model);

    //constructed graph with edges all the s_i set to true,
    // iterate over all its vertices
    // and for each pair of edges, check whether there is also the third edge in the triangle -> if not, conflict
    // runtime n max_d^2, where max degree should be relatively small, since we look only at all s_ij that are true
    for (int vertex_j = 0; vertex_j < num_vertices; vertex_j++) {
        //next, iterate over pairs of neighbours
        for(auto i_it = pairing_graph[vertex_j].begin(); i_it != pairing_graph[vertex_j].end(); ++i_it){
            for (auto k_it = std::next(i_it); k_it != pairing_graph[vertex_j].end(); ++k_it) {
                //from the pairing graph we know that s_ij and s_jk are true, check for s_ik
                if(not get_model_ij(model, *i_it, *k_it)){
                    //found a conflict: s_ij and s_jk are true but s_ik is not
                    conflicts.emplace_back(*i_it, vertex_j, *k_it);
                }
            }
        }
    }
    if(conflicts.empty()){
        assert(find_conflicts_naive().empty());
    }
    return conflicts;
}


ConflictList IncSatGC::find_conflicts_paper() const {
    //store found conflicts
    ConflictList conflicts;
    Model model = get_model();
    int num_colors;
    switch (options.strategy) {
        case Options::TopDown:
            num_colors = upper_bound;
            break;
        case Options::BottomUp:
            num_colors = lower_bound;
            break;
        case Options::SingleK:
            assert(options.specific_num_colors.has_value());
            num_colors = options.specific_num_colors.value();
            break;
        default:
            if(options.verbosity >= Options::Debug){
                std::cout << "c Debug: Number of colors used for checker was determined from model\n";
            }
            num_colors = get_num_colors_from_model();
            break;
    }

    // c[] data structure that stores set of available colors for vertex v in c[v].
    // Encoded as c[v][i] = 1 if color i is available; all colors are available in the beginning
    std::vector<Bitset> c(num_vertices, Bitset(num_colors).set());
    //main loop of check algorithm
    for (int u = 0; u < num_vertices; ++u) {
        auto index = c[u].find_first();
        if (index == Bitset::npos) {//no available color found
            //this case is fine, can be that c[u] is empty because it was involved in an earlier conflict
            continue;
        }
        //set c[u] = {index}, i.e. choose color for u
        c[u].reset();
        c[u].set(index);
        //inner loop
        for (int v = u + 1; v < num_vertices; ++v) {
            //check whether -s_uv is in the pairing
            if (not get_model_ij(model, u, v)) {
                //update c[v] = c[v]\c[u] i.e. set c[v][index] = false
                c[v].reset(index);
                if (c[v].none()) {
                    for (int w = 0; w < u; ++w) {
                        if( get_model_ij(model, w, v) and (c[w] == c[u])) {
                            //CONFLICT TYPE ii, add transitivity(u,w,v)
                            conflicts.emplace_back(u, w, v); //w < u < v
                        }
                    }
                    //it's fine if we cannot find any 'w' here
                }
            } else { //s_uv is in the pairing
                //update c[v] = c[v] intersect c[u]
                c[v] &= c[u];
                if (c[v].none()) {
                    bool w_exists = false; //just to make sure such a 'w' always exists
                    for (int w = 0; w < u; ++w) {
                        if(c[w].none()){
                            //do not consider if w is empty, otherwise redundant triples are added
                            continue;
                        }
                        bool model_wv = get_model_ij(model, w, v);
                        if (model_wv and (c[w] != c[u])) {
                            //CONFLICT TYPE iii, add transitivity(u,v,w)
                            conflicts.emplace_back(u, v, w);
                            w_exists = true;
                        }
                        else if( not model_wv and (c[w] == c[u])) {
                            //CONFLICT TYPE i, add transitivity(w,u,v)
                            conflicts.emplace_back(w, u, v);
                            w_exists = true;
                        }
                    }
                    assert(w_exists);
                }
            }
        }
    }
    if(conflicts.empty()){
        assert(find_conflicts_naive().empty());
    }
    return conflicts;
}


void IncSatGC::add_conflict_clauses(const ConflictList &conflicts) {
    for (const auto &conflict : conflicts) {
        auto [i,j,k] = conflict;
        add_transitivity(i,j,k);
    }
}


int IncSatGC::do_zykov_propagator() {
    if (options.solver != Options::CaDiCaL) {
        throw std::runtime_error("Propagators are only supported by Cadical");
    }
    switch (options.strategy) {
        case Options::SingleK:
            zykov_propagator_single_k();
            break;
        case Options::TopDown:
            zykov_propagator_top_down();
            break;
        case Options::BottomUp:
            zykov_propagator_bottom_up();
            break;
        default:
            throw std::runtime_error("Invalid option for search strategy.");
    }
    return 0;
}

void IncSatGC::init_zykov_propagator() {
    //new cadical propagator stuff
    // 1. Connect externalPropagator to cadical solver
    std::shared_ptr<CaDiCaLAdaptor::Solver> cadical_extended = std::dynamic_pointer_cast<CaDiCaLAdaptor::Solver>(solver);
    assert(cadical_extended != nullptr);
    CaDiCaL::Solver* cadical_solver = &(cadical_extended->solver);
    assert(cadical_solver != nullptr);
    zykov_propagator = std::make_unique<CadicalZykovPropagator>(CadicalZykovPropagator(*this));
    cadical_solver->connect_external_propagator(zykov_propagator.get());


    // 2. add every edge variable s_ij as observed variable to the propagator
    assert(not complement_graph_adjacency.empty());
    for (int i = 0; i < num_vertices; ++i) {
        for (int j : complement_graph_adjacency[i] ) {
            if(i > j){
                //only want to look at (i,j) once
                continue;
            }
            //add var s_ij to the observed vars of the propagator
            // std::cout << "Add s_" << i << "," << j << " (" << sij_indices.get(i,j) + 1 << ") as observed var\n";
            cadical_solver->add_observed_var(sij_indices.get(i,j) + 1);//index + 1 for cadical
        }
    }
    // also add c_j vars as observed
    for (int j = 1; j < num_vertices; ++j) {
        // std::cout << "Add c_" << j << " (" << c_indices[j - 1] + 1 << ") as observed var\n";
        cadical_solver->add_observed_var(c_indices[j - 1] + 1);//index + 1 for cadical
    }
    //3. do not add any transitivity constraints, these will be handled by the propagator
    // the cardinality constraints also follow externally, or later potentially as part of the propagator too
}

bool IncSatGC::zykov_propagator_single_k() {
    return zykov_encoding_single_k();
}

int IncSatGC::zykov_propagator_top_down() {
    return zykov_encoding_top_down();
}

int IncSatGC::zykov_propagator_bottom_up() {
    return zykov_encoding_bottom_up();
}


int IncSatGC::original_paper_configuration() {
    //first step: build zykov encoding, including unit clauses!
    stats.start_phase(Statistics::BuildEncoding);
    //initialise indices of s_ij and c_j in member fields
    //we have variables s_i,j for all i,j in V, only consider i<j so that is n(n-1)/2 in upper triangular matrix
    //only look at i,j such that (i,j) is not an edge
    //and variables c_j for 1 <= j <= n-1. keep track of indices in vectors
    sij_indices = UpperTriangle(num_vertices);
    c_indices.resize(num_vertices - 1);
    int var_index = 0;

    complement_graph_adjacency = graph.get_complement().get_neighbor_list();
    for (int i = 0; i < num_vertices; ++i) {
        for (int j = i + 1; j < num_vertices; ++j) {
            //i and j are not adjacent in original graph
            sij_indices.set(i,j, var_index);
            var_index++;
        }
    }
    stats.num_sij_vars = var_index;

    for (int j = 1; j < num_vertices; ++j) {
        c_indices[j - 1] = var_index;
        var_index++;
    }
    stats.num_cj_vars = num_vertices - 1;
    assert(stats.num_cj_vars == var_index - stats.num_sij_vars);

    add_vars(c_indices.back() + 1);//since we start counting at 0 we have 1 more variable than last index

    //here we do add unit literals for all the edges
    Graph::NeighborList neighbors = graph.get_neighbor_list();
    for (int i = 0; i < num_vertices; ++i) {
        for (int j : neighbors[i]) {
            if(i > j) {
                continue;
            }
            solver->addClause(mkLit(sij_indices.get(i, j), true));
        }
    }

    //next are the definitions of c_j as in equation 2
    assert(not tmp_clause.size());
    vec<Lit> second_clause;
    for (int j = 1; j < num_vertices ; ++j) {
        //for each j, add clause c_j <-> AND -s_ij
        //which in cnf is AND(-c_j v -s_ij) and (ORs_ij v c_j)
        //do first clauses while collecting literals for second
        second_clause.capacity(j);
        for (int i = 0; i < j; ++i) {
            if(sij_indices.get(i, j) == -1){
                //is an edge, thus leave it out
                continue;
            }
            tmp_clause.push(mkLit(sij_indices.get(i, j), true));
            tmp_clause.push(mkLit(c_indices[j - 1], true));
            solver->addClause(tmp_clause);
            tmp_clause.clear();

            second_clause.push(mkLit(sij_indices.get(i, j), false));
        }
        second_clause.push(mkLit(c_indices[j - 1], false));
        solver->addClause(second_clause);
        second_clause.clear();
    }


    stats.end_phase(Statistics::BuildEncoding);

    int num_colors = lower_bound; //1toN in paper, actually starts at 2
    //also set lower_bound to 1 and assert that we have no upper bound
    assert(lower_bound == 2 and upper_bound == graph.ncount());

    //done with building encoding, start bottom up loop
    while(lower_bound != upper_bound) {//outer loop, increasing num colors
        add_incremental_at_most_k(num_colors - 1); //only encoding that works, even if not the one from paper
        print_sat_size();
        bool res = search_add_conflicts_and_solve(); //inner loop, cegar approach of adding new clauses until no more conflicts
        notify_new_bound(res, num_colors);
        if (res){//found a k-coloring, became satisfiable
            notify_upper_bound(num_colors);
            current_best_coloring = obtain_coloring_from_model();
            assert(*std::max_element(current_best_coloring.begin(), current_best_coloring.end()) + 1 == num_colors);
            assert(lower_bound == upper_bound);
            //loop stops here
        }
        else {//not satisfiable, increase number of colors and continue
            notify_lower_bound(num_colors + 1);
            num_colors++;
        }
    }
    //chromatic number is lower_bound=upper_bound
    stats.solved = true;
    return upper_bound;
}


void IncSatGC::register_write_cleanup_on_exit() const{
    auto exitcode = [](int sign) {
            std::cout << "Aborting due to signal " << sign << "\n";
            write_and_cleanup_on_exit();
            std::cout.flush();
            std::exit(128 + sign);
        };

        //^C
        signal(SIGINT, exitcode);
        // abort()
        signal(SIGABRT, exitcode);
        // sent by "kill" command
        signal(SIGTERM, exitcode);
        //^Z
        signal(SIGTSTP, exitcode);
        // ulimit -t
        signal(SIGXCPU, exitcode);
        // ulimit -v
        signal(SIGSEGV, exitcode);
}

void IncSatGC::write_and_cleanup() {
    stats.end_phase(Statistics::Total);
    if(options.verbosity >= Options::Normal){
        stats.print_stats();
    }
    else if(options.verbosity == Options::Quiet and (lower_bound == upper_bound)){
        std::cout << "o "<< lower_bound << "\n";
    }
    if(not options.stats_csvfile.empty()){
        stats.write_stats();
    }
    if(not options.coloringfilepath.empty()){
        write_optimal_coloring();
    }
}






