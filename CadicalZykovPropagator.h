#ifndef INCSATGC_CADICALPROPAGATOR_H
#define INCSATGC_CADICALPROPAGATOR_H


#include <iostream>
#include <memory>
#include <vector>
#include <queue>
#include <map>

#include <boost/dynamic_bitset.hpp>
#include "cadical.hpp"

#include "GraphMatrix.h"
#include "Statistics.h"

using ConflictQueue = std::deque< std::vector<int> >;
using Bitset = boost::dynamic_bitset<>;

//declare IncSatGC
struct UpperTriangle;
class IncSatGC;

class CadicalZykovPropagator : public CaDiCaL::ExternalPropagator{

public:
    bool is_lazy = false; // lazy propagator only checks complete assignments
    bool are_reasons_forgettable = true; // Reason external clauses can be deleted
    bool is_tainting = true; // The external clauses must trigger restore (unless frozen)

    ~CadicalZykovPropagator () override ;

    // Notify the propagator about assignments to observed variables.
    // The notification is not necessarily eager. It usually happens before
    // the call of propagator callbacks and when a driving clause is leading
    // to an assignment.
    //
    // void notify_assignment (int lit, bool is_fixed) override ;
    void notify_assignment (const std::vector<int>& lits) override;
    void notify_new_decision_level () override ;
    void notify_backtrack (size_t new_level) override ;

    // Check by the external propagator the found complete solution (after
    // solution reconstruction). If it returns false, the propagator must
    // provide an external clause during the next callback.
    //
    bool cb_check_found_model (const std::vector<int> &model) override ;

    // Ask the external propagator for the next decision literal. If it
    // returns 0, the solver makes its own choice.
    //
    int cb_decide () override ;

    // Ask the external propagator if there is an external propagation to make
    // under the current assignment. It returns either a literal to be
    // propagated or 0, indicating that there is no external propagation under
    // the current assignment.
    //
    int cb_propagate () override ;

    // Ask the external propagator for the reason clause of a previous
    // external propagation step (done by cb_propagate). The clause must be
    // added literal-by-literal closed with a 0. Further, the clause must
    // contain the propagated literal.
    //
    int cb_add_reason_clause_lit (int propagated_lit) override ;

    // The following two functions are used to add external clauses to the
    // solver during the CDCL loop. The external clause is added
    // literal-by-literal and learned by the solver as an irredundant
    // (original) input clause. The clause can be arbitrary, but if it is
    // root-satisfied or tautology, the solver will ignore it without learning
    // it. Root-falsified literals are eagerly removed from the clause.
    // Falsified clauses trigger conflict analysis, propagating clauses
    // trigger propagation. In case chrono is 0, the solver backtracks to
    // propagate the new literal on the right decision level, otherwise it
    // potentially will be an out-of-order assignment on the current level.
    // Unit clauses always (unless root-satisfied, see above) trigger
    // backtracking (independently from the value of the chrono option and
    // independently from being falsified or satisfied or unassigned) to level
    // 0. Empty clause (or root falsified clause, see above) makes the problem
    // unsat and stops the search immediately. A literal 0 must close the
    // clause.
    //
    // The external propagator indicates that there is a clause to add.
    // The parameter of the function allows the user to indicate that how
    // 'forgettable' is the external clause. Further, there is a member-parameter
    // called 'is_tainting' that indicates if the external clauses are irredundant
    // (almost always they are, see below the rare exceptions).
    // There are overall four possible scenarios defined by these two Boolean flags:
    // - Non-Forgettable & Tainting Clause [default value].
    //    - Policy to Forget: The clause will not be deleted, unless it is
    //      directly implied by some of the other irredundant clauses (e.g. root
    //      satisfied).
    //    - Tainting: The clause will be considered as irredundant and thus the
    //      literals of it will be tainted and might trigger restore steps
    //      (except frozen literals).
    //    - Example: It is always safe and correct to use this clause type.
    // - Forgettable & Tainting Clause.
    //    - Policy to Forget: The clause will be considered for deletion during
    //      clause database reduction rounds.
    //    - Tainting: The clause will be considered as irredundant and thus the
    //      literals of it will be tainted and might trigger restore steps (except
    //      frozen literals).
    //    - Example usage: Theory lemmas in SMT solvers. Anything that is
    //      irredundant from the SAT perspective (hence must be considered for
    //      tainting), but allowed to be deleted because the external propagator
    //      (e.g. theory solver) can always re-derive it.
    // - Forgettable & Not-Tainting Clause.
    //    - Policy to Forget: The clause will be considered for deletion during
    //      clause database reduction rounds.
    //    - Tainting: The clause is assumed to be derivable in propositional
    //      logic from the current set of irredundant and redundant clauses (i.e.
    //      (\phi \land \rho) \implies C, as in the side condition of Learn- in
    //      incremental inprocessing), therefore the literals of it are not
    //      tainted and do not trigger restore steps.
    //    - Example usage: Use it ONLY to share learned clauses between parallel
    //      SAT solvers. It is NOT correct to use this clause type unless the
    //      clause is indeed directly derivable by the SAT solver alone.
    // - Non-Forgettable & Non-Tainting Clause:
    //    - Policy to Forget: The clause will not be deleted, unless it is
    //      directly implied by some of the other irredundant clauses (e.g. root
    //      satisfied).
    //    - Redundancy: The clause is assumed to be derivable in propositional
    //      logic from the current set of irredundant and redundant clauses (i.e.
    //      (\phi \land \rho) \implies C, as in the side condition of Learn- in
    //      incremental inprocessing), therefore the literals of it are not
    //      tainted and do not trigger restore steps.
    //    - Example usage: Same as previous case, the only difference is that the
    //      clauses will be kept around (can be useful in cases where their
    //      derivation was expensive).
    //
    //
    // Reason clauses of external propagation steps are assumed to be forgettable,
    // and tainting by default. Use parameter 'reason_forgettable' to change it.
    //
    // DISCLAIMER: Any tainting clause is correct to be added ONLY if it is
    // guaranteed that the SAT solver uses only RUP-based clause addition
    // techniques (so no BVA) and the literals are either tainted or frozen.
    // See our Incremental Inprocessing SAT'19 paper for more details on it.
    //
    // The naming of the clause types will might change later (to coordinate
    // with IPASIR-2).
    //
    bool cb_has_external_clause (bool& is_forgettable) override;

    // The actual function called to add the external clause.
    //
    int cb_add_external_clause_lit () override ;


    /*
     * New section for new data members and functions
     */

    //need reference to IncSatGC instance to access data and call its functions, in particular the find conflict ones
    IncSatGC* INSTANCE; //use raw pointer to not take ownership of the superclass
    //constructor that constructs reference to INSTANCE
    explicit CadicalZykovPropagator(IncSatGC& reference);

    Statistics& stats;
    const Options& options;

    //some helpful member fields which are mostly the ones of the INSTANCE
    int num_vertices;
    int num_colors;
    void update_num_colors(int num_colors_);
    UpperTriangle* sij_indices;
    [[nodiscard]] inline int getCadicalIndex(int i, int j) const;
    //start with simple functionality to check whether the returned model is correct
    // and otherwise return conflict clauses when asked later on
    // for now, use Triangle Checker for finding conflicts
    [[nodiscard]] ConflictQueue find_conflicts_in_model(const std::vector<int> &model) const;

    //for the transitivity propagation i need two big things:
    // a copy of the graph that I modify, i.e. add edges/contract vertices
    // and a union find kinda structure that keeps tracks of the bags and updates them when to vertices are merged
    MGraph mgraph;



    int highest_sij_var;
    int highest_var;
    int num_assigned = 0;
    int current_level = 0;
    //store all found propagations and corresponding reason clauses
    std::deque<int> propagations;
    std::deque< std::vector<int> > clauses;
    //once a literal has been propagated, remember the reason clause in literal_to_clause, either positive or negative
    std::vector<std::vector<int>> literal_to_clause_pos;
    std::vector<std::vector<int>> literal_to_clause_neg;
    //store the decisions and propagations done, i.e. the trail
    std::vector< std::vector<int> > current_trail;
    //store the assignment of variables so far, done with a signed char* of vals and negatively indexed as in cadical
    size_t vsize = 0; // actually allocated variable data size
    int max_var = 0;  // internal maximum variable index
    signed char *vals = nullptr; // assignment [-max_var,max_var]
    void enlarge_vals(int new_max_var);
    [[nodiscard]] signed char val(int lit) const;
    void set_val(int lit, signed char val);
    static int sign (const int lit) { return (lit > 0) - (lit < 0); } //branchless sign computation
    //vector in which to store any external clauses we want to add
    std::vector<std::vector<int>> external_clauses;


    //main functions called when notified of new decision
    void propagate_same(int u_old, int v_old);
    void propagate_differ(int u_old, int v_old);
    //bitsets as members to avoid new allocation every time
    Bitset nu_without_nv;
    Bitset nv_without_nu;
    //helper functions to propagate the decisions
    // u,v are merged and v,w are merged, propagate that uw are merged and give transitivity(u,v,w) as a reason
    void propagate_positive(int u, int v, int w);
    // u,v are merged and v,w are adjacent!, propagate that uw are adjacent and give transitivity(v,u,w) as a reason
    void propagate_negative(int u, int v, int w);

    //functions for decision strategy of next literal
    int next_decision();
    [[nodiscard]] int first_literal() const;
    [[nodiscard]] int ISUN_literal() const;
    [[nodiscard]] std::tuple<int, int> get_edge_candidates(const Bitset &clique) const;
    int imitate_dsatur_literal();
    [[nodiscard]] int bag_size_literal() const;

    //functions and fields to compute and store maximal cliques, often not of maximum size
    std::vector<Bitset> maximal_cliques;
    int max_clique_size;
    void compute_cliques();
    [[nodiscard]] bool need_to_recompute_cliques() const;
    //some helper functions to determine whether cliques need to be recomputed and which bounds we compute
    long long num_assignments_last_clique_computation = 0; //need to recompute cliques if assignments were done
    [[nodiscard]] bool compute_clique_clauses() const;
    [[nodiscard]] bool compute_mycielsky_clauses() const;
    bool first_call_after_backtrack; //only call expensive mycielsky bound after a backtrack
    //for cliques to not mess up the bottom-up approach,
    //we need an activation literal to disable clauses for cliques of smaller sizer later on
    int bottom_up_clique_assumption_variable = 0;
    [[nodiscard]] bool need_bottom_up_clique_assumption_variable() const;
    [[nodiscard]] bool bottom_up_clique_assumption_variable_is_set() const;

    //helper function to put clique explanation into a vector
    [[nodiscard]] std::vector<int> clique_explanation_clause(const Bitset &clique) const;
    //function that checks cliques for being of size > num_colors, and adds their explanation
    void check_for_clique_clauses();
    //given a clique, produce the explanation clause and add it to external_clauses
    void add_clique_explanation_clause(const Bitset &clique);
    //function that tries to compute mycielsky extension based off of cliques
    void check_for_mycielsky_clauses();
    //given a mycielsky bound, produce the explanation clause and add it to external_clauses
    void add_mycielsky_explanation_clause(const MGraph::SubGraph &subgraph);


    //use external binary CliSAT to find a strong clique on the current subgraph
    // modifies max_clique_size and stores clique in max_cliques
    void find_cliques_externally();

    //functions for clique-based pruning as pre-processing of the subproblem graph
    // and to compute dominated vertices as a "good" decision to be made
    //store vertices whose neighbourhood changed, only those are relevant for neighbourhood search
    std::vector<int> touched_vertices;
    //only a decision, return the lit that was decided on
    int find_dominated_vertex_decisions();
    void find_vertex_fusions();
    void find_edge_additions();
    void find_clique_based_pruning();

    //functions to look for a coloring on graph of current node in Zykov tree
    [[nodiscard]] bool compute_coloring() const;
    void check_for_coloring() const;

};


#endif //INCSATGC_CADICALPROPAGATOR_H
