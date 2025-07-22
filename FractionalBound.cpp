#include "FractionalBound.h"

double fractional_chromatic_number_exactcolors(const std::vector<Bitset>& adj){
    fflush(nullptr);  // flush all stdio buffers
    int saved_out  = dup(STDOUT_FILENO);
    int saved_err  = dup(STDERR_FILENO);
    int nullfd     = open("/dev/null", O_WRONLY|O_CLOEXEC);
    dup2(nullfd, STDOUT_FILENO);
    dup2(nullfd, STDERR_FILENO);

    // convert graph from vector of bitsets to n,m and array of edges
    const int n = static_cast<int>(adj.size());
    // flatten upper-triangle into (u,v) edge list
    std::vector<int> edges;
    edges.reserve(n * (n - 1) / 2);
    for (int u = 0; u < n; ++u) {
        for (int v = u + 1; v < n; ++v) {
            if (adj[u].test(v)) {
                edges.push_back(u);
                edges.push_back(v);
            }
        }
    }
    const int m = static_cast<int>(edges.size() / 2);

    // build ExactColors problem
    COLORproblem node_problem;
    COLORproblem_init(&node_problem);
    COLORproblem_init_with_graph(&node_problem, n, m, edges.data());

    COLORparms* p = &node_problem.parms;
    COLORparms_init(p);
    COLORset_dbg_lvl(0);

    int frac = -1.0;
    if (!compute_lower_bound(&node_problem.root_cd, &node_problem)) {
        frac = node_problem.root_cd.lower_bound;
    }

    //reset stdout
    dup2(saved_out, STDOUT_FILENO);
    dup2(saved_err, STDERR_FILENO);
    close(saved_out); close(saved_err); close(nullfd);

    //Clean up
    COLORproblem_free(&node_problem);
    // COLORlp_free_env();

    return frac;
}