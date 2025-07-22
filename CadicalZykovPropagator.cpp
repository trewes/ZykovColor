#include "CadicalZykovPropagator.h"

#include "IncSatGC.h" //make type information available to use here for pointer of IncSatGC object

# define  DEBUG_PRINT(X) //do {X} while(0)  //debug print that can be disabled easily

# define  PROP_TIMING(X) //do {X} while(0)  //call to timing functions that can easily be disabled

CadicalZykovPropagator::~CadicalZykovPropagator() = default;

void CadicalZykovPropagator::notify_assignment(const std::vector<int>& lits) {
    PROP_TIMING(stats.start_phase(Statistics::PropagatorAssignment););
    stats.prop_num_assignments += static_cast<int>(lits.size());
    DEBUG_PRINT(std::cout << "Was notified of assignments " << lits << " at level " << current_level << "\n";);
    for(int lit : lits) {
        int absLit = abs(lit);
        if(val(lit) != 0){
            DEBUG_PRINT(std::cout << "propagated lit " << lit << " already has an assignment " << val(lit) << "\n";);
            //might get notified about the same assignment twice
            assert(val(lit) == (lit > 0 ? 1 : -1));
            continue;
        }
        set_val(lit, sign (lit));
        current_trail.back().push_back(lit);
        num_assigned++;

        if(absLit > highest_sij_var) {
            DEBUG_PRINT(
                if(absLit > highest_var){//higher than encoding vars, must be one of the (possibly old) assumption literals
                    assert(need_bottom_up_clique_assumption_variable());
                    std::cout << "Was notified of assignment of activation literal " << lit << " at level " << current_level << "\n";
                }
                else
                    std::cout << "Was notified of assignment of literal " << lit << " = c_" << (std::abs(lit)  - INSTANCE->c_indices.front())
                          << " at level " << (current_level) << "\n";
                );
            //do not propagate anything for assignments of cj literals or activation literal
            continue;
        }

        auto [u,v] = sij_indices->get_ij(std::abs(lit) - 1);
        int u_rep = mgraph.vertex_rep[u];
        int v_rep = mgraph.vertex_rep[v];
        DEBUG_PRINT(
            std::cout << "Was notified of assignment of literal " << lit << " = {" << u << "," << v << "} represented by "
                      << u_rep << ", " << v_rep << " at level " << current_level << "\n";
        );

        if(mgraph.has_edge(u_rep,v_rep) or mgraph.is_contracted(u_rep,v_rep)){
            //if assignment was already made on the graph side, ignore it here
            //this also ignores contradicting assignments but they are still handled in the propagations
            continue;
        }
        //handle same (s_u,v = true) and different (s_u,v = false) operations on graph
        if(lit > 0){
            propagate_same(u,v);
        }
        else {
            propagate_differ(u,v);
        }
    }
    PROP_TIMING(stats.end_phase(Statistics::PropagatorAssignment););
}

void CadicalZykovPropagator::notify_new_decision_level() {
    assert(current_level + 1 == current_trail.size());
    current_level++;
    DEBUG_PRINT(std::cout << "Was notified of new decision level " << current_level << "\n";);
    //normally there should be no more propagations but sometimes full trail is renotified
    //which also increases decision level without asking propagations (but they are assigned so its fine)
    propagations.clear();
    clauses.clear();
    assert(external_clauses.empty());
    maximal_cliques.clear();
    max_clique_size = -1;

    current_trail.emplace_back();
    mgraph.notify_new_level();
    assert(mgraph.added_edges_limiter.size() == mgraph.current_level);
    stats.prop_max_level = std::max(stats.prop_max_level, current_level);
    touched_vertices.clear();
	//pre-emptive resize of stats vector for larger levels
	if(current_level >= stats.prop_node_depth_history.size()){
		int new_size = current_level * 2;
		stats.prop_node_depth_history.resize(new_size, 0);
		if(options.use_clique_explanation_clauses){
			stats.prop_clique_pruning_level.resize(new_size,0);
		}
		if(options.use_mycielsky_explanation_clauses){
			stats.prop_myc_pruning_level.resize(new_size,0);
		}
		if(options.enable_positive_pruning){
			stats.prop_positive_pruning_level.resize(new_size,0);
		}
		if(options.enable_negative_pruning){
			stats.prop_negative_pruning_level.resize(new_size,0);
		}
	}
	stats.prop_node_depth_history[current_level]++;
}

void CadicalZykovPropagator::notify_backtrack(size_t new_level) {
    DEBUG_PRINT(std::cout << "Was notified of backtrack: " << new_level << "\n";);
    PROP_TIMING(stats.start_phase(Statistics::PropagatorBacktrack););
    stats.prop_num_backtracks++;
    int length = current_level - new_level;
    if(length >= stats.prop_backtrack_size.size()) {
        stats.prop_backtrack_size.resize(2 * length, 0);
    }
    stats.prop_backtrack_size[length]++;
    if(options.enable_detailed_backtracking_stats) { //collect current level, level we backtrack to, and reason
        stats.prop_detailed_backtrack_list.insert(stats.prop_detailed_backtrack_list.end(), {current_level, new_level, stats.backtrack_reson});
        stats.backtrack_reson = 0;
    }

    assert(current_level > new_level);
    while (current_level > new_level){
        for(auto l_it = current_trail.back().rbegin(); l_it != current_trail.back().rend(); ++l_it){
            set_val(*l_it, 0);
            num_assigned--;
        }
        current_trail.pop_back();
        current_level--;
        assert(current_level + 1 == current_trail.size());
    }

    mgraph.notify_backtrack_level(static_cast<int>(new_level));
    assert(mgraph.added_edges_limiter.size() == mgraph.current_level);
    DEBUG_PRINT(mgraph.print(););

    //reset all propagation info
    propagations.clear();
    clauses.clear();
    //also reset clique explanations
    external_clauses.clear();
    maximal_cliques.clear();
    max_clique_size = -1;
    first_call_after_backtrack = true;
    touched_vertices.clear();
    PROP_TIMING(stats.end_phase(Statistics::PropagatorBacktrack););
}

bool CadicalZykovPropagator::cb_check_found_model(const std::vector<int> &model) {
    DEBUG_PRINT(std::cout << "EP checked found model\n";);
    assert(clauses.empty());
    clauses = find_conflicts_in_model(model);
    if(not clauses.empty()){
        std::cout << "cb_check_found_model : false, because checker found conflicts\n";
        throw std::runtime_error("Checked model should be correct?");
    }
    //model was correct which means a valid coloring was found, compare to time when heuristic found one, if any
    if(stats.heuristic_found_coloring) {
        stats.heuristic_theoretical_time_improvement += stats.current_total_time() - stats.heuristic_coloring_time_point;
        stats.heuristic_found_coloring = false;
    }
    return true;
}

int CadicalZykovPropagator::cb_decide() {
    DEBUG_PRINT(std::cout << "EP was asked for next decision literal\n";);
    PROP_TIMING(stats.start_phase(Statistics::PropagatorDecide););
    int decision_literal = 0;
    if(options.use_dominated_vertex_decisions) {
        //returns the literal to contract two dominated vertices or 0 if none was found
        decision_literal = find_dominated_vertex_decisions();
    }
    if(decision_literal == 0) {
        decision_literal = next_decision();
    }
    stats.prop_num_decisions++;
    PROP_TIMING(stats.end_phase(Statistics::PropagatorDecide););
    return decision_literal;
}

int CadicalZykovPropagator::cb_propagate() {
    DEBUG_PRINT(std::cout << "EP was asked if there is an external propagation\n";);
    assert(clauses.size() == propagations.size());
    if(propagations.empty()) {
        //done with transitivity propagations for now, check if there are prunings to be propagated
        find_clique_based_pruning();
    }
    while (not propagations.empty()){
        //get next propagation and store corresponding reason clause in literal_to_clause
        int literal = propagations.front();
        assert(literal != 0);
        int abslit = std::abs(literal);
        propagations.pop_front();
        auto reason_clause = clauses.front();
        clauses.pop_front(); // delete clause
        if(val(literal) == sign(literal)) {
            return cb_propagate();//literal was already assigned to value of propagation, take next propagation
        }
        //make sure the literal is part of the given clause
        assert(std::find(reason_clause.begin(), reason_clause.end(), literal) != reason_clause.end());
        if(literal > 0) {
            literal_to_clause_pos[abslit] = reason_clause;
        }
        else {
            literal_to_clause_neg[abslit] = reason_clause;
        }

        DEBUG_PRINT(
            int i = std::get<0>(sij_indices->get_ij(std::abs(literal) - 1));
            int j = std::get<1>(sij_indices->get_ij(std::abs(literal) - 1));
            std::cout << "cb_propagate " << literal << " {" << i << "," << j << "}\n";
        );
        stats.prop_num_propagations++;
        return literal;
    }
    assert(propagations.empty() and clauses.empty());
    return 0;
}

int CadicalZykovPropagator::cb_add_reason_clause_lit(int propagated_lit) {
    //given a previously propagated literal, produce the reason literal by literal now
    DEBUG_PRINT(
        int i = std::get<0>(sij_indices->get_ij(std::abs(propagated_lit) - 1));
        int j = std::get<1>(sij_indices->get_ij(std::abs(propagated_lit) - 1));
        std::cout << "EP was asked for the next lit reason of external propagation " << propagated_lit
                  << "{" << i << "," << j << "} \n";
    );
    int abslit = std::abs(propagated_lit);
    assert(propagated_lit != 0);
    auto &reason = (propagated_lit > 0 ? literal_to_clause_pos[abslit] : literal_to_clause_neg[abslit]);
    if (not reason.empty()){
        auto lit = reason.back();
        reason.pop_back();
        DEBUG_PRINT(std::cout << "cb_add_reason_clause_lit " << lit << "\n";);
        return lit;
    }
    //gave reason clause, increment stat
    stats.prop_num_reason_clauses++;
    return 0;
}

bool CadicalZykovPropagator::cb_has_external_clause(bool& is_forgettable) {
    DEBUG_PRINT(std::cout << "EP was asked if there are external clauses\n";);
    //only look for external clauses if there are no more propagations to be made to keep assignemnts and graph in sync
    if(not propagations.empty() or not bottom_up_clique_assumption_variable_is_set()) {
        return false;
    }
    if(compute_clique_clauses()) {
        check_for_clique_clauses();
    }
    if(compute_mycielsky_clauses()) {
        check_for_mycielsky_clauses();
    }
    if(compute_coloring()){
        check_for_coloring();
    }
    if (compute_fractional_bound()) {
        check_for_fractional_bound();
    }
    first_call_after_backtrack = false;

    if(not external_clauses.empty()){
        DEBUG_PRINT(std::cout << "cb_has_external_clause : true, because external_clauses not empty\n";);
        is_forgettable = true; //the explanation clauses can be forgotten, we only really want to use them to backtrack once
        return true;
    }
    return false;
}

int CadicalZykovPropagator::cb_add_external_clause_lit() {
    DEBUG_PRINT(std::cout << "EP was asked for the next external clause literal\n";);
    std::vector<int> &clause = external_clauses.back();
    if (not clause.empty()){
        int lit = clause.back();
        clause.pop_back();
        DEBUG_PRINT(std::cout << "cb_add_external_clause_lit " << lit << "\n";);
        return lit;
    }
    // delete last clause
    external_clauses.pop_back();
    stats.prop_num_external_clauses++;
    return 0;
}


/*
 * New section for new data members and functions
 */

CadicalZykovPropagator::CadicalZykovPropagator(IncSatGC &reference)
    : stats(reference.stats), options(reference.options)
{
    INSTANCE = &reference;
    num_vertices = INSTANCE->num_vertices;
    num_colors = num_vertices;
    sij_indices = &INSTANCE->sij_indices;
    mgraph = MGraph(num_vertices, INSTANCE->graph.ecount(), INSTANCE->graph.elist());

    highest_sij_var = INSTANCE->c_indices.front(); //no +1 so this is the last sij var
    highest_var = INSTANCE->c_indices.back() + 1; //+1 because of cadical indices
    propagations = {};
    clauses = {};
    literal_to_clause_pos = std::vector<std::vector<int>>(highest_sij_var + 1);
    literal_to_clause_neg = std::vector<std::vector<int>>(highest_sij_var + 1);
    current_trail.emplace_back();
    enlarge_vals(highest_var + 1);//since lits and thus indices are positive, not starting at zero

    external_clauses = {};
    maximal_cliques = {};
    max_clique_size = {};
    first_call_after_backtrack = false;

	//pre-allocate space for some of the stat vectors

    stats.prop_node_depth_history.resize(500,0);
    stats.prop_backtrack_size.resize(50,0);
    if(options.enable_detailed_backtracking_stats) {
        stats.prop_detailed_backtrack_list.reserve(100000);
    }
    if(options.use_clique_explanation_clauses){
    	stats.prop_clique_pruning_level.resize(500,0);
    }
    if(options.use_mycielsky_explanation_clauses){
	    stats.prop_myc_pruning_level.resize(500,0);
    }
    if(options.enable_positive_pruning){
    	stats.prop_positive_pruning_level.resize(500,0);
    }
    if(options.enable_negative_pruning){
    	stats.prop_negative_pruning_level.resize(500,0);
    }

    flag_fractional_timed_out = INSTANCE->flag_fractional_timed_out;
}

void CadicalZykovPropagator::update_num_colors(int num_colors_) {
    num_colors = num_colors_;
    //extra considerations when using bottom up strategy with clique explanation clauses
    //have to add clique clauses with activation literal and assume literal to be false,
    //then add literal true as unit clause to deactivate the clique clauses when increasing num_colors
    if(need_bottom_up_clique_assumption_variable()) {
        auto solver = std::dynamic_pointer_cast<CaDiCaLAdaptor::Solver>(INSTANCE->solver);
        //assumption literal exists, add unit clause
        if(bottom_up_clique_assumption_variable != 0) {
            assert(bottom_up_clique_assumption_variable > 0);
            solver->solver.add(bottom_up_clique_assumption_variable);
            solver->solver.add(0);
        }
        //find first variable index we can use as assumption variables
        if(bottom_up_clique_assumption_variable == 0) {
            //one variable after last, can use it as new variable for assumption
            bottom_up_clique_assumption_variable = solver->solver.vars() + 1;
        }
        else {
            //previously chose an assumption variable, can just increment by one
            bottom_up_clique_assumption_variable++;
        }
        //add variable as observed so that it can be used in propagation reason clauses
        solver->solver.add_observed_var(bottom_up_clique_assumption_variable);
        //have to check whether assumption var is assigned and safe to use, i.e. that enough space for all variables is there
        if(bottom_up_clique_assumption_variable > max_var) {
            //should never be necessary but just to be safe
            enlarge_vals(bottom_up_clique_assumption_variable);
        }
        assert(bottom_up_clique_assumption_variable <= vsize);
        //assume it as false so rest of clique clause has to be satisfied
        solver->solver.assume(-bottom_up_clique_assumption_variable);
    }
}

int CadicalZykovPropagator::getCadicalIndex(int i, int j) const {
    //cadical numbering starts at 1, other variable numbering starts at 0 though
    return sij_indices->get(i,j) + 1;
}

ConflictQueue CadicalZykovPropagator::find_conflicts_in_model(const std::vector<int> &model) const {
    ConflictQueue conflicts;
    Graph::NeighborList pairing_graph(num_vertices);
    //add edge for each s_ij set to true in model, iterate over each non-edge to check
    for (int i = 0; i < num_vertices; ++i) {
        for (int j : INSTANCE->complement_graph_adjacency[i] ) {
            if(i > j){
                //only want to look at (i,j) once
                continue;
            }
            if (model[sij_indices->get(i,j)] > 0) {
                pairing_graph[i].push_back(j);
                pairing_graph[j].push_back(i);
            }
        }
    }
    //constructed graph with edges all the s_i set to true, iterate over all its vertices
    // and for each pair of edges, check whether there is also the third edge in the triangle -> if not, conflict
    // runtime n max_d^2, where max degree should be relatively small, since we look only at all s_ij that are true
    for (int vertex_j = 0; vertex_j < num_vertices; vertex_j++) {
        //next, iterate over pairs of neighbours
        for (auto vertex_i = pairing_graph[vertex_j].begin(); vertex_i != pairing_graph[vertex_j].end(); ++vertex_i) {
            for (auto vertex_k = std::next(vertex_i); vertex_k != pairing_graph[vertex_j].end(); ++vertex_k) {
                bool model_ik = (sij_indices->get(*vertex_i,*vertex_k) == -1) ? false : (model[sij_indices->get(*vertex_i,*vertex_k)]  > 0);
                if(not model_ik){
                    //found a conflict: s_ij and s_jk are true but s_ik is not
                    std::vector<int> conflict;
                    conflict.push_back( - (getCadicalIndex(*vertex_i, vertex_j)));
                    conflict.push_back( - (getCadicalIndex(vertex_j, *vertex_k)));
                    if(sij_indices->get(*vertex_i,*vertex_k) != -1){
                        conflict.push_back(   getCadicalIndex(*vertex_i, *vertex_k));
                    }
                    DEBUG_PRINT(std::cout << "Conflict vertices" << *vertex_i << " " << vertex_j << " " << *vertex_k << ", literals " << conflict.back() << "\n";);
                    conflicts.push_back(conflict);
                }
            }
        }
    }
    return conflicts;
}

void CadicalZykovPropagator::enlarge_vals(int new_max_var) {
    if(new_max_var < vsize) {
        //no need to make vals[] larger, there is enough space
        assert(new_max_var > max_var);
        max_var = new_max_var;
        return;
    }
    size_t new_vsize = 2 * new_max_var;
    signed char *new_vals;
    const size_t bytes = 2u * new_vsize;
    new_vals = new signed char[bytes]; // g++-4.8 does not like ... { 0 };
    memset (new_vals, 0, bytes);
    new_vals += new_vsize; //shift so we can access at negative indices

    if (vals) {
        memcpy (new_vals - max_var, vals - max_var, 2u * max_var + 1u);
        vals -= vsize;
        delete[] vals;
    }
    vals = new_vals;
    vsize = new_vsize;
    max_var = new_max_var;
}

signed char CadicalZykovPropagator::val(int lit) const {
    assert (-max_var <= lit);
    assert (lit);
    assert (lit <= max_var);
    return vals[lit];
}

void CadicalZykovPropagator::set_val(int lit, signed char val) {
    assert (-1 <= val);
    assert (val <= 1);
    assert (-max_var <= lit);
    assert (lit);
    assert (lit <= max_var);
    vals[lit] = val;
    vals[-lit] = val; //cadical assigns -val and then takes care when retreiving the value, we want the assignment regardless of sign!
}

void CadicalZykovPropagator::propagate_same(int u_old, int v_old) {
    assert(not mgraph.has_edge(u_old, v_old));
    assert(not mgraph.is_contracted(u_old, v_old));
    DEBUG_PRINT(std::cout << "propagate_same " << u_old << " " << v_old << "\n";);

    //get representative vertex, but also keep originals for correct propagation
    int u_rep = mgraph.vertex_rep[u_old];
    int v_rep = mgraph.vertex_rep[v_old];
    assert(not mgraph.has_edge(u_rep, v_rep));
    assert(not mgraph.is_contracted(u_rep, v_rep));

    //get different neighbours of u and v for adding edges later
    nu_without_nv = MGraph::setminus(mgraph.gmatrix[u_rep], mgraph.gmatrix[v_rep]);
    nv_without_nu = MGraph::setminus(mgraph.gmatrix[v_rep], mgraph.gmatrix[u_rep]);

    //vertex u and v are merged, so are vertex u and v' for all v' in in bag[v]
    // all v and u' for u' in bag[u], and finally, u'v' after that
    for(int vp : mgraph.bag[v_rep]){
        if(vp == v_old) { continue;}
        propagate_positive(u_old, v_old, vp);
    }
    for(int up : mgraph.bag[u_rep]){
        if(up == u_old) { continue;}
        propagate_positive(v_old, u_old, up);
    }
    for(int vp : mgraph.bag[v_rep]){
        if(vp == v_old) { continue;}
        for(int up : mgraph.bag[u_rep]){
            if(up == u_old) { continue;}
            propagate_positive(vp, u_old, up);
            propagate_positive(vp, v_old, up);
        }
    }
    //for neighbours w of u but not of v, separate v and w
    for (int vp : mgraph.bag[v_rep]) {
        for (int w =  nu_without_nv.find_first(); w != Bitset::npos; w = nu_without_nv.find_next(w)) {
            propagate_negative(vp, u_old, w);
        }
    }
    //likewise for neighbours w of v but not of u, separate u and w
    for (int up : mgraph.bag[u_rep]) {
        for (int w =  nv_without_nu.find_first(); w != Bitset::npos; w = nv_without_nu.find_next(w)) {
            propagate_negative(up, v_old, w);
        }
    }

    mgraph.contract_vertices(u_rep, v_rep);
    //keep track of vertices whose neighbourhood changed, in this case the vertex which was merged into
    assert(mgraph.vertex_rep[std::min(u_rep, v_rep)] == std::min(u_rep, v_rep)); //merged into smaller of the two
    touched_vertices.push_back(std::min(u_rep, v_rep));
}

void CadicalZykovPropagator::propagate_differ(int u_old, int v_old) {
    assert(not mgraph.has_edge(u_old, v_old));

    DEBUG_PRINT(std::cout << "propagate_differ " << u_old << " " << v_old << "\n";);

    //get representative vertex, but also keep originals for correct propagation
    int u_rep = mgraph.vertex_rep[u_old];
    int v_rep = mgraph.vertex_rep[v_old];
    if(mgraph.has_edge(u_rep,v_rep)){
        return; //there is an edge, already separated
    }

    //vertex u and v are separated, so are vertex u and v' for all v' in in bag[v]
    // all v and u' for u' in bag[u], and finally, u'v' after that
    for(int up : mgraph.bag[u_rep]){
        if(up == u_old) { continue;}
        propagate_negative(up, u_old, v_old);
    }
    for(int vp : mgraph.bag[v_rep]){
        if(vp == v_old) { continue;}
        propagate_negative(vp, v_old, u_old);
    }
    for(int up : mgraph.bag[u_rep]){
        if(up == u_old) { continue;}
        for(int vp : mgraph.bag[v_rep]){
            if(vp == v_old) { continue;}
            propagate_negative(up, u_old, vp);
            // propagate_negative(vp, v_old, up); //i dont need this
        }
    }

    //add edge between representatives and bags in vertices
    mgraph.separate_vertices(u_rep, v_rep);
    //keep track of vertices whose neighbourhood changed, in this case both vertices gained an edge
    touched_vertices.push_back(u_rep);
    touched_vertices.push_back(v_rep);
}

void CadicalZykovPropagator::propagate_positive(int u, int v, int w) {
    assert(u != w && v != w);
    assert(u != v); //should already have checked this before calling the function
    assert( sij_indices->get(u,w) != -1 );
    assert(not mgraph.is_contracted(u,w));
    assert(not mgraph.has_edge(u,w));

    int propagated_lit = getCadicalIndex(u, w);
    DEBUG_PRINT(std::cout << "propagate_positive " << u << " " << v << " " << w << " i.e. " << propagated_lit << " = {" << u << "," << w << "}\n";);
    assert(val(propagated_lit) == 0);

    //want to propagate s_u,w = true with reason -s_uv v -s_vw v s_uw
    propagations.push_back(propagated_lit);
    std::vector<int> new_clause;
    new_clause.push_back(-getCadicalIndex(u, v));
    new_clause.push_back(-getCadicalIndex(v, w));
    new_clause.push_back(propagated_lit);
    clauses.push_back(new_clause);

    DEBUG_PRINT(
        std::cout << "Added " << propagated_lit << " = {" << u << "," << w << "}" << " to propagation queue with reason ";
        for (auto lit : new_clause){
            std::cout << lit << " ";
        }std::cout <<"\n";
    );
}

// u and v are merged, v,w are adjacent, propagate edge u-w
void CadicalZykovPropagator::propagate_negative(int u, int v, int w) {
    assert(u != w && v != w);
    assert(u != v); //should already have checked this before calling the function
    assert(sij_indices->get(u,w) != -1);//original graph already has an edge
    assert(not mgraph.has_edge(u,w));

    int propagated_lit = -getCadicalIndex(u, w); // propagate edge between u and w
    DEBUG_PRINT(std::cout << "propagate_negative " << u << " " << v << " " << w << " i.e. " << propagated_lit << " = {" << u << "," << w << "}\n";);
    assert(val(propagated_lit) == 0);

    //want to propagate s_u,w = false with reason -s_uv v s_vw v -s_uw
    propagations.push_back(propagated_lit);
    std::vector<int> new_clause;
    new_clause.push_back(-getCadicalIndex(u, v));
    if(sij_indices->get(v,w) != -1){ //no edge in original graph is what we want to check here
        new_clause.push_back(getCadicalIndex(v, w));
    }
    new_clause.push_back(propagated_lit);
    clauses.push_back(new_clause);

    DEBUG_PRINT(
        std::cout << "Added " << propagated_lit << " = {" << u << "," << w << "}" << " to propagation queue with reason ";
        for (auto lit : new_clause){
            std::cout << lit << " ";
        }std::cout <<"\n";
    );
}


int CadicalZykovPropagator::next_decision() {
    switch (options.zykov_propagator_decision_strategy) {
        case Options::CadicalZykov:
            return 0; //let cadical decide
        case Options::FirstLiteral:
            return first_literal();
        case Options::ISUN:
            return ISUN_literal();
        case Options::ImitateDsatur:
            return imitate_dsatur_literal();
        case Options::BagSize:
            return bag_size_literal();
        default:
            throw std::runtime_error("Invalid decision strategy or Zykov propagator.");
    }
}

int CadicalZykovPropagator::first_literal() const {
    //simply return the first unset edge variable
    for (auto iIt = mgraph.vertices.begin(); iIt != mgraph.vertices.end(); ++iIt){
        for (auto jIt = std::next(iIt); jIt != mgraph.vertices.end(); ++jIt){
            int i = *iIt;
            int j = *jIt;
            if(mgraph.has_edge(i,j)) {
                continue;
            }
            assert(not mgraph.is_contracted(i,j)); //should not be contracted if both are still in vertex set
            int decision_literal = getCadicalIndex(i,j);
            DEBUG_PRINT(std::cout << "FirstLiteral decided on " << decision_literal << " = {" << i << "," << j << "}\n";);
            return decision_literal; //what about negative literal assignment?
        }
    }
    //possibly find no unset edge variable but there are still cardinality constraint variables to assign, let cadical do this
    return 0;
}

int CadicalZykovPropagator::ISUN_literal() const {
    //ISUN strategy: {u,v} not in E s.t. d(u)+d(v) is maximal
    int choice_u = -1;
    int choice_v = -1;
    int max_sum = 0;
    //only iterate over vertices whose neighbourhood changed, only their degree increased
    for(int i : (touched_vertices.empty() ? mgraph.vertices : touched_vertices)){
        if(not mgraph.nodeset[i]) {
            //vertex might have been contracted into another
            continue;
        }
        for (int j : INSTANCE->complement_graph_adjacency[i]){ //can only dominate each other if not adjacent
            if((not mgraph.nodeset[j]) or mgraph.has_edge(i,j)) {
                continue;
            }
            assert(mgraph.vertex_rep[j] == j and i != j);
            assert(not mgraph.is_contracted(i,j)); //should not be contracted if both are still in vertex set
            int ij_sum = static_cast<int>(mgraph.gmatrix[i].count() + mgraph.gmatrix[j].count());
            if( ij_sum > max_sum) {
                choice_u = i;
                choice_v = j;
                max_sum = ij_sum;
            }
        }
    }
    //found an unset variable and set max_sum, as well as choice u,v
    if (max_sum > 0) {
        assert(choice_u >= 0 and choice_v >= 0);
        int decision_literal = getCadicalIndex(choice_u,choice_v);
        DEBUG_PRINT(std::cout << "ISUN decided on " << decision_literal << " = {" << choice_u << "," << choice_v << "}\n";);
        return decision_literal; //what about negative literal assignment?
    }
    //possibly find no unset edge variable but there are still cardinality constraint variables to assign, let cadical do this
    return 0;
}


std::tuple<int, int> CadicalZykovPropagator::get_edge_candidates(const Bitset &clique) const {
    //from clique, find branching candidates
    // i.e. v adjacent to most vertices in clique and u in clique not adjacent to v
    int max_colors = -1;
    int max_degree = -1;
    int max_v = -1;
    //choose v from remaining vertices not in clique
    if(mgraph.nodeset == clique) {
        //catch case that graph is only the clique and there are no other vertices left
        return {-1, -1};
    }
    for (int v : mgraph.vertices){
        if(clique[v]) {
            continue;
        }
        //compute size of intersection of neighbours of v and clique
        int colors = (mgraph.gmatrix[v] & clique).count();
        assert(colors < clique.size());
        if(colors < max_colors) {
            continue;
        }
        //also include degree as criteria
        int degree = mgraph.gmatrix[v].count();

        if(colors > max_colors or degree > max_degree ) {
            //found vertex that improves number of colors, or has the same colors but larger degree
            max_colors = colors;
            max_degree = degree;
            max_v = v;
        }
    }
    assert(max_v >= 0); //check that we found a vertex
    //found a v adjacent to most clique vertices (and largest degree), next look for an adjacent u
    //u has to be in clique and not a neighbour of max_v
    int max_u = (clique - mgraph.gmatrix[max_v]).find_first();
    if(max_u == Bitset::npos) {
        //catch case that graph is not connected, i.e. v is not adjacent to clique and no u exists
        return {-1, -1};
    }
    assert(max_u >= 0);
    assert(not mgraph.has_edge(max_u, max_v));


    DEBUG_PRINT(std::cout << "Chosen u,v: " << max_u << " and " << max_v << " on level " << current_level << "\n";);

    //candidate will be max_u, max_v
    return {max_u, max_v};
}

int CadicalZykovPropagator::imitate_dsatur_literal() {
    //decision is based on cliques, so we need to make sure we computed them
    if (need_to_recompute_cliques()){
        compute_cliques();
    }
    //go through cliques and find a literal to return as decision
    for(const Bitset &clique : maximal_cliques) {
        auto [choice_u,choice_v] = get_edge_candidates(clique);
        if(choice_u == -1 and choice_v == -1) {
            //no vertices were found from clique choice, try next clique
            continue;
        }
        //valid choice was found, return decision literal
        assert(choice_u >= 0 and choice_v >= 0);
        int decision_literal = getCadicalIndex(choice_u,choice_v);
        DEBUG_PRINT(std::cout << "ImitateDsatur decided on " << decision_literal << " = {" << choice_u << "," << choice_v << "}\n";);
        return decision_literal;
    }
    //if no clique was successful, let cadical decide
    DEBUG_PRINT(std::cout << "ImitateDsatur decided on " << 0 << "\n";);
    return 0;
}

int CadicalZykovPropagator::bag_size_literal() const {
    //choose vertices u,v such that the sum of their bag sizes is maximal
    int choice_u = -1;
    int choice_v = -1;
    int max_sum = 0;
    int degree_sum = 0;
    //only iterate over vertices whose neighbourhood changed
    for(int i : (touched_vertices.empty() ? mgraph.vertices : touched_vertices)){
        if(not mgraph.nodeset[i]) {
            //vertex might have been contracted into another
            continue;
        }
        for (int j : INSTANCE->complement_graph_adjacency[i]){ //can only dominate each other if not adjacent
            if((not mgraph.nodeset[j]) or mgraph.has_edge(i,j)) {
                continue;
            }
            assert(mgraph.vertex_rep[j] == j and i != j);
            assert(not mgraph.is_contracted(i,j)); //should not be contracted if both are still in vertex set
            int ij_sum = static_cast<int>(mgraph.bag[i].size() + mgraph.bag[j].size());
            if( max_sum < ij_sum) {
                choice_u = i;
                choice_v = j;
                max_sum = ij_sum;
                degree_sum = static_cast<int>(mgraph.gmatrix[i].count() + mgraph.gmatrix[j].count());
            }
            else if(ij_sum == max_sum) {
                if(degree_sum < static_cast<int>(mgraph.gmatrix[i].count() + mgraph.gmatrix[j].count())) {
                    choice_u = i;
                    choice_v = j;
                    max_sum = ij_sum;
                    degree_sum = static_cast<int>(mgraph.gmatrix[i].count() + mgraph.gmatrix[j].count());
                }
            }
        }
    }
    //found an unset variable and set max_sum, as well as choice u,v, ony use it if bag sizes are non-trivial
    if (max_sum > 2) {
        assert(choice_u >= 0 and choice_v >= 0);
        int decision_literal = getCadicalIndex(choice_u,choice_v);
        DEBUG_PRINT(std::cout << "BagSize decided on " << decision_literal << " = {" << choice_u << "," << choice_v << "}\n";);
        return decision_literal; //what about negative literal assignment?
    }

    //possibly find no unset edge variable but there are still cardinality constraint variables to assign, let cadical do this
    return 0;
}

void CadicalZykovPropagator::compute_cliques() {
    PROP_TIMING(stats.start_phase(Statistics::PropagatorComputeCliques););
    maximal_cliques.clear();
    max_clique_size = -1;
    num_assignments_last_clique_computation = stats.prop_num_assignments;
    //run any clique algorithm to store largest cliques in maximal_cliques
    max_clique_size = mgraph.greedy_cliques(maximal_cliques, options.prop_clique_limit);
    if (options.mnts_length > 0 and max_clique_size <= num_colors){
        std::vector<Bitset> mnts_cliques;
        int mnts_size = mgraph.mnts_clique(mnts_cliques, num_colors + 1,
                                options.mnts_length, options.mnts_depth,
                                static_cast<int>(stats.prop_num_clique_computations)); //use changing but deterministic random seed
        if (mnts_size > max_clique_size) {
            max_clique_size = mnts_size;
            maximal_cliques = mnts_cliques;
        }
    }
    assert(not maximal_cliques.empty());
    stats.prop_num_clique_computations++;
    stats.prop_num_maximal_cliques_computed += static_cast<int>(maximal_cliques.size());
    if(max_clique_size == num_colors) {
        stats.prop_num_tight_cliques_computed += static_cast<int>(maximal_cliques.size());
    }
    PROP_TIMING(stats.end_phase(Statistics::PropagatorComputeCliques););
}


bool CadicalZykovPropagator::need_to_recompute_cliques() const {
    //cliques might not yet exist or be valid anymore if new assignments were made on the same level
    return maximal_cliques.empty() or (num_assignments_last_clique_computation != stats.prop_num_assignments);
}

bool CadicalZykovPropagator::compute_clique_clauses() const {
    assert(propagations.empty());
    //options is set and we are done with assumptions, actual assignments have started
    return options.use_clique_explanation_clauses and num_assigned > INSTANCE->encoder_assumptions.size();
}

bool CadicalZykovPropagator::compute_mycielsky_clauses() const {
    assert(propagations.empty());
    //options is set, we have backtracked before this call and we are done with assumptions, actual assignments have started
    return options.use_mycielsky_explanation_clauses and first_call_after_backtrack and external_clauses.empty()
            and num_assigned > INSTANCE->encoder_assumptions.size();
}

bool CadicalZykovPropagator::need_bottom_up_clique_assumption_variable() const {
    return options.strategy == Options::BottomUp
            and (options.use_clique_explanation_clauses
                or options.enable_positive_pruning
                or options.enable_negative_pruning);
}

bool CadicalZykovPropagator::bottom_up_clique_assumption_variable_is_set() const {
    //either don't care about activation literal or ceck that it is assigned
    assert(options.strategy != Options::BottomUp or bottom_up_clique_assumption_variable > 0);
    assert(options.strategy != Options::BottomUp or bottom_up_clique_assumption_variable <= max_var);
    //can potentially be 1 too after fixing it via unit clause, then just return that it is not set (to what we need)
    return options.strategy != Options::BottomUp or val(bottom_up_clique_assumption_variable) == -1;
}

std::vector<int> CadicalZykovPropagator::clique_explanation_clause(const Bitset &clique) const {
    std::vector<int> clause;
    if(need_bottom_up_clique_assumption_variable()) {
        assert(bottom_up_clique_assumption_variable > 0);
        assert(bottom_up_clique_assumption_variable_is_set());
        //add assumed literal to clique explanation so that the clause can be deactivated later
        clause.push_back(bottom_up_clique_assumption_variable);
    }
    //for u,v in clique, add OR e_r(u),r(v)
    for (int u = clique.find_first(); u != Bitset::npos; u = clique.find_next(u)) {
        for (int v = clique.find_next(u); v != Bitset::npos; v = clique.find_next(v)) {
            //clique vertices should be represented by right vertex
            assert(mgraph.vertex_rep[u] == u);
            assert(mgraph.vertex_rep[v] == v);
            assert(mgraph.has_edge(u,v));
            //find lit for u,v and add e_u,v to external clause list
            int elit = getCadicalIndex(u,v);
            if (elit == 0) {
                //u,v was an edge in original graph and cannot be contracted, leave it out of clause
                continue;
            }
            assert(val(elit) == -1);
            clause.push_back(elit);
        }
    }
    return clause;
}

void CadicalZykovPropagator::check_for_clique_clauses() {
    if (not bottom_up_clique_assumption_variable_is_set()) {
        return;
    }
    if(need_to_recompute_cliques()) {
        compute_cliques();
    }
    assert(not maximal_cliques.empty());
    PROP_TIMING(stats.start_phase(Statistics::PropagatorCliqueClauses););
    //check for cliques of size > num_colors and add clauses that one of the edges of the clique has to be contracted to avoid a new color
    //(because a clique of that size would need more colors than we are currently allowing)
    if(max_clique_size > num_colors) {
        for(const auto & clique : maximal_cliques) {
            //clique is too large and one edge should have been contracted instead, add the reason clause for this
            add_clique_explanation_clause(clique);
			stats.prop_clique_pruning_level[current_level]++;
        }
    }
    PROP_TIMING(stats.end_phase(Statistics::PropagatorCliqueClauses););
}


void CadicalZykovPropagator::add_clique_explanation_clause(const Bitset &clique) {
    DEBUG_PRINT(std::cout << "Trying to explain clique " << clique << " of size " << clique.count() << "\n";);
    assert(mgraph.is_clique(clique));
    external_clauses.push_back(clique_explanation_clause(clique));
    stats.prop_num_clique_successes++;
    stats.backtrack_reson = 1;
    DEBUG_PRINT(std::cout << "Gave explanation " << external_clauses.back() << " for clique " << clique << " of size " << clique.count() << "\n";);
}

void CadicalZykovPropagator::check_for_mycielsky_clauses() {
    if (not bottom_up_clique_assumption_variable_is_set()) {
        return;
    }
    if(need_to_recompute_cliques()) {
        compute_cliques();
    }
    assert(not maximal_cliques.empty());
    PROP_TIMING(stats.start_phase(Statistics::PropagatorMycielskyClauses););
    assert(max_clique_size <= num_colors); //otherwise clique pruning would have been done instead
    int gap = num_colors - max_clique_size;
    if(gap < options.mycielsky_threshold) {
        for(const auto & clique : maximal_cliques) {
            if(gap >= stats.mycielsky_calls.size()) {
                stats.mycielsky_calls.resize(gap + 1, 0);
            }
            stats.mycielsky_calls[gap]++;
            //start algorithm with clique subgraph; given as bitset and build graph from that
            MGraph::SubGraph subgraph(clique);
            //try to extend subgraph to next generalised mycielsky, need at least gap + 1 succesful iterations
            int iterations = mgraph.mycielsky_extension_clique(subgraph, gap + 1);
            if(iterations > gap) { //only explain subgraph if the proven bound was good enough to cross threshold
                assert(iterations == gap + 1);
                //graph is too large and one edge should have been contracted instead, add the reason clause for this
                add_mycielsky_explanation_clause(subgraph);
                if(gap >= stats.mycielsky_sucesses.size()) {
                    stats.mycielsky_sucesses.resize(gap + 1, 0);
                }
                stats.mycielsky_sucesses[gap]++;
				stats.prop_myc_pruning_level[current_level]++;
                stats.backtrack_reson = 2;
            }
        }
    }
    PROP_TIMING(stats.end_phase(Statistics::PropagatorMycielskyClauses););
}

void CadicalZykovPropagator::add_mycielsky_explanation_clause(const MGraph::SubGraph &subgraph) {
    DEBUG_PRINT(std::cout << "Trying to explain mycielsky subgraph " << subgraph.nodes << " of size " << subgraph.num_vertices << "\n";);
    external_clauses.emplace_back();
    if(need_bottom_up_clique_assumption_variable()) {
        assert(bottom_up_clique_assumption_variable > 0);
        assert(bottom_up_clique_assumption_variable_is_set());
        //add assumed literal to clique explanation so that the clause can be deactivated later
        external_clauses.back().push_back(bottom_up_clique_assumption_variable);
    }
    //for u,v in subgraph, add OR e_r(u),r(v)
    for (int u : subgraph.nodes) {
        for (int v = subgraph.matrix[u].find_next(u); v != Bitset::npos; v = subgraph.matrix[u].find_next(v)) {
            //mycielsky subgraph vertices should be represented by right vertex
            assert(mgraph.vertex_rep[u] == u);
            assert(mgraph.vertex_rep[v] == v);
            assert(mgraph.has_edge(u,v));
            //find lit for u,v and add -e_u,v to external clause list
            int elit = getCadicalIndex(u,v);
            if (elit == 0) {
                //u,v was an edge in original graph and cannot be contracted, leave it out of clause
                continue;
            }
            assert(val(elit) == -1); //might not be true if propagation has not happened yet!
            external_clauses.back().push_back(elit);
        }
    }
    DEBUG_PRINT(std::cout << "Gave explanation " << external_clauses.back() << " for mycielsky subgraph " << subgraph.nodes << " of size " << subgraph.num_vertices << "\n";);
}


//helper function to produce subgraph description string in dimacs format
std::string get_subgraph_dimacs(const MGraph &subgraph) {
    std::stringstream ss;
    int edge_count = 0;

    for (size_t i = 0; i < subgraph.vertices.size(); ++i) {
        int u = subgraph.vertices[i];
        for (size_t j = i + 1; j < subgraph.vertices.size(); ++j) {
            int v = subgraph.vertices[j];
            if (subgraph.has_edge(u,v)) {
                ss << " e " << u+1 << " " << v+1 << "\n";
                edge_count++;
            }
        }
    }
    std::stringstream header;
    header << "p edge " << subgraph.size << " " << edge_count << "\n";
    return header.str() + ss.str();
}


int CadicalZykovPropagator::find_dominated_vertex_decisions() {
    //touched vertices contains all vertices who had their neighbourhood changed, either through contraction or adding an edge
    for (int i : touched_vertices){
        if(not mgraph.nodeset[i]) {
            //vertex might have been contracted into another
            continue;
        }
        for (int j : INSTANCE->complement_graph_adjacency[i]){ //can only dominate each other if not adjacent
            if((not mgraph.nodeset[j]) or mgraph.has_edge(i,j)) {
                continue;
            }
            assert(mgraph.vertex_rep[j] == j and i != j);
            //vertex i cannot become dominated by edge additions, so only need to check whether it dominates j
            if(mgraph.gmatrix[j].is_subset_of(mgraph.gmatrix[i])) {
                //j is dominated by i since N(i) contains all of N(j) or vice versa
                int lit = getCadicalIndex(i,j);
                assert(val(lit) == 0);
                DEBUG_PRINT(
                std::cout << "DomVertex decided on " << lit  << " since "<< i << " and " << j << " dominate each other"
                          << " on level " << current_level << "\n";
                );
                stats.prop_num_dominated_vertex_decisions++;
                return lit;
            }
        }
    }
    return 0;
}


void CadicalZykovPropagator::find_vertex_fusions() {
    PROP_TIMING(stats.start_phase(Statistics::PropagatorPositivePruning););
    assert(max_clique_size == num_colors);
    //test for vertex fusion, i.e., u in C and non-adjacent v in V s.t. v is adjacent to all of C\u, then we can merge u,v
    for(auto & clique : maximal_cliques) {
        bool exit = false;
        for (int u = clique.find_first(); u != Bitset::npos and not exit; u = clique.find_next(u)) {
            clique.reset(u); //for easier adjacency check with vertices v, intersection has to be clique minus u
            for (auto v : INSTANCE->complement_graph_adjacency[u]) {
                if(not mgraph.nodeset[v] or mgraph.has_edge(u,v)) {
                    continue;
                }
                if(clique.is_subset_of(mgraph.gmatrix[v])) {
                    //in this case, we can merge u and v
                    int propagated_lit = getCadicalIndex(u, v);
                    propagations.push_back(propagated_lit);
                    //as reason clause, we give all edges of the clique and edges from v to the clique
                    clique.set(u); //include u in clique explanation!
                    std::vector<int> clause = clique_explanation_clause(clique);
                    clique.reset(u);
                    for (int w = clique.find_first(); w != Bitset::npos; w = clique.find_next(w)) {
                        assert(mgraph.vertex_rep[w] == w);
                        assert(mgraph.has_edge(v,w));
                        //find lit for u,v and add e_u,v to external clause list
                        int elit = getCadicalIndex(v,w);
                        if (elit == 0) {
                            //u,v was an edge in original graph and cannot be contracted, leave it out of clause
                            continue;
                        }
                        assert(val(elit) == -1);
                        clause.push_back(elit);
                    }
                    clause.push_back(propagated_lit);
                    clauses.push_back(clause);
                    stats.prop_positive_prunings++;
					stats.prop_positive_pruning_level[current_level]++;
                    DEBUG_PRINT(std::cout << "Vertex fusion propagates " << propagations.back() << " merge " << u << ","
                        << v << " with reason " << clauses.back() << " on level " << current_level
                        << " clique " << clique << " and N(v) " << mgraph.gmatrix[v] << "\n";);
                    exit = true;//only want to find one conflict per clique, as a balance with speed
                    break;
                }
            }
            clique.set(u);
        }
    }
    PROP_TIMING(stats.end_phase(Statistics::PropagatorPositivePruning););
}

void CadicalZykovPropagator::find_edge_additions() {
    PROP_TIMING(stats.start_phase(Statistics::PropagatorNegativePruning););
    assert(max_clique_size == num_colors);
    //test for possible edge addition, i.e., for u,v in V\C, test that for all w, (u,w) or (v,w) in E
    for(int u : touched_vertices){
        if(not mgraph.nodeset[u]) {
            continue;
        }
        // bool exit = false;
        for (auto v : INSTANCE->complement_graph_adjacency[u]) {
            if(u > v or not mgraph.nodeset[v] or mgraph.has_edge(u,v)) {
                continue;
            }
            for(auto & clique : maximal_cliques) {
                if(clique[u] or clique[v]) {
                    continue;
                }
                //got u,v in V\C, now check clique for condition C subset (N(u) U N(v))
                if(clique.is_subset_of(mgraph.gmatrix[u] | mgraph.gmatrix[v])) {
                    //in this case, we can add an edge between u and v
                    int propagated_lit = -getCadicalIndex(u, v);
                    propagations.push_back(propagated_lit);
                    //as reason clause, we give all edges of the clique and edges from x and y to the clique
                    std::vector<int> clause = clique_explanation_clause(clique);
                    for (int w = clique.find_first(); w != Bitset::npos; w = clique.find_next(w)) {
                        int vertex = (mgraph.has_edge(u,w) ? u : v);
                        assert((vertex == u and mgraph.has_edge(u,w)) or (vertex == v and mgraph.has_edge(v,w)));
                        assert(mgraph.vertex_rep[w] == w);
                        assert(mgraph.has_edge(vertex,w));
                        //find lit for u,v and add e_u,v to external clause list
                        int elit = getCadicalIndex(vertex,w);
                        if (elit == 0) {
                            //vertex,w was an edge in original graph and cannot be contracted, leave it out of clause
                            continue;
                        }
                        assert(val(elit) == -1);
                        clause.push_back(elit);
                    }
                    clause.push_back(propagated_lit);
                    clauses.push_back(clause);
                    stats.prop_negative_prunings++;
					stats.prop_negative_pruning_level[current_level]++;
                    DEBUG_PRINT(std::cout << "Edge addition propagates " << propagations.back() << " edge " << u << ","
                            << v << " with reason " << clauses.back() << " on level " << current_level
                            << " clique " << clique << " and N(u) " << mgraph.gmatrix[u] << " and N(v) " << mgraph.gmatrix[v] << "\n";);
                    // exit = true;
                    break;
                }
            }
            // if(exit){break;} //continue with next vertex u when a conflict for u,v to limit time spent here
        }
    }
    PROP_TIMING(stats.end_phase(Statistics::PropagatorNegativePruning););
}

void CadicalZykovPropagator::find_clique_based_pruning() {
    if(not bottom_up_clique_assumption_variable_is_set()) {
        //assumption not set yet, can't use it in bottom-up propagations so skip these for now
        return;
    }
    if(need_to_recompute_cliques() and (options.enable_positive_pruning or options.enable_negative_pruning)) {
        compute_cliques();
    }
    //prunings only hold if clique is same size as num_colors we are looking for
    if(max_clique_size != num_colors) {
        return;
    }
    //check for positive and negative pruning, i.e., vertex fusions and edge additions
    if(options.enable_positive_pruning) {
        find_vertex_fusions();
    }
    if(options.enable_negative_pruning and first_call_after_backtrack and propagations.empty()) {
        find_edge_additions();
    }
    //since these are prunings that need to hold (otherwise problem is not k-colorable),
    //it is okay to propagate all of them at the same time, even if they conflict each other.
    //with correct handling, this simply leads to a conflict and backtrack
}

bool CadicalZykovPropagator::compute_coloring() const {
    assert(propagations.empty());
    //options is set, we have backtracked before this call and we are done with assumptions, actual assignments have started
    return options.zykov_coloring_algorithm != Options::None and first_call_after_backtrack and external_clauses.empty()
            and num_assigned > INSTANCE->encoder_assumptions.size();
}

void CadicalZykovPropagator::check_for_coloring() const {
    //try to compute better coloring on changed graph, only returns size and coloring is not used
    PROP_TIMING(stats.start_phase(Statistics::TestColorTime););
    int coloring_size = stats.upper_bound;
    switch (options.zykov_coloring_algorithm) {
        case Options::FastDsatur:
            coloring_size = mgraph.dsatur_coloring((maximal_cliques.empty() ? Bitset(num_vertices) : maximal_cliques.front()));
            break;
        case Options::SortedSEQ:
            coloring_size = mgraph.sequential_coloring((maximal_cliques.empty() ? Bitset(num_vertices) : maximal_cliques.front()));
            break;
        case Options::IteratedIS:
            coloring_size = mgraph.IS_extract((maximal_cliques.empty() ? Bitset(num_vertices) : maximal_cliques.front()));
            break;
        case Options::IteratedSEQ:
            coloring_size = mgraph.ISEQ();
            break;
        case Options::None:
            throw std::runtime_error("Should have previously checked that no coloring algorithm was to be used.");
    }

    if(coloring_size <= num_colors) {
        std::cout << "Found an improved coloring using " << coloring_size << " colors!!!!! at " << INSTANCE->stats.current_total_time()
        << " density of " << mgraph.density() <<"\n";
        if(not stats.heuristic_found_coloring) { //only update first time coloring is found
            stats.heuristic_found_coloring = true;
            stats.heuristic_coloring_time_point = stats.current_total_time();
        }
    }
    PROP_TIMING(stats.end_phase(Statistics::TestColorTime););
}


bool CadicalZykovPropagator::compute_fractional_bound() const {
    if (flag_fractional_timed_out) {
        //computation for initial bound computation took too long, do not use bounding in propagator
        return false;
    }
    if (options.use_fractional_bound and first_call_after_backtrack and external_clauses.empty()
            and num_assigned > INSTANCE->encoder_assumptions.size()) {
        if (mgraph.density() >= options.frac_density) { //density threshold criteria
            return true;
        }
    }
    return false;
}

void CadicalZykovPropagator::check_for_fractional_bound() {
    PROP_TIMING(stats.start_phase(Statistics::PropagatorFractionalBound););
    if (options.verbosity >= Options::Debug) {
        std::cout << "Checking fractional bound on level " << current_level << " with "
                  << std::accumulate(current_trail.begin(), current_trail.end(), 0, [](int sum, auto& v){return sum + v.size();}) << " assignments and "
                  << std::accumulate(current_trail.begin(), current_trail.end(), 0, [&](int sum, auto& v){return sum +
                      std::count_if(v.begin(), v.end(), [&](int a){return a > 0 and a - 1 < highest_sij_var;});}) << " merged"
                  << " at density " << mgraph.density() << "\n";
    }
    stats.fractional_bound_calls++;
    auto start = Statistics::cpuTime();
    //build graph on which to compute bound
    std::vector<int> index_mapping;
    std::vector<Bitset> active_subgraph = mgraph.get_active_subgraph(index_mapping);
    //compute bound
    double frac = fractional_chromatic_number_exactcolors(active_subgraph);
    auto total = Duration(Statistics::cpuTime() - start);
    stats.full_fractional_time += total;

    //if bound causes a conflict, add an external clause
    if (frac > num_colors) {
        stats.fractional_bound_success++;
        if (options.verbosity >= Options::Debug) {
            std::cout << std::setprecision(25) << "PRUNE! Fractional bound: " << frac << " " << total << "\n";
        }
        std::vector<int> tmp;
        // add trivial reason clause to cause backtrack
        assert(current_level + 1 == current_trail.size());
        for (int i = 1; i <= current_level; i++) {
            int decision_lit = current_trail[i][0];
            if (decision_lit != bottom_up_clique_assumption_variable and decision_lit - 1 > highest_sij_var) {
                continue;
            }
            assert(std::dynamic_pointer_cast<CaDiCaLAdaptor::Solver>(INSTANCE->solver)->solver.is_decision(decision_lit));
            tmp.push_back(-decision_lit);
        }
        external_clauses.push_back(tmp);
    }
    else {
        if (options.verbosity >= Options::Debug) {
            std::cout << std::setprecision(25) << "Fractional bound: " << frac << " " << total << "\n";
        }
    }
    PROP_TIMING(stats.end_phase(Statistics::PropagatorFractionalBound););
}




