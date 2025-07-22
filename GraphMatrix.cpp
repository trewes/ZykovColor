#include "GraphMatrix.h"


MGraph::MGraph(const int n, const int m, const std::vector<int> &edge_list) {
    assert(2*m == edge_list.size());
    size = n;
    num_vertices = n;
    gmatrix = std::vector<Bitset>(size, Bitset(size));
    nodeset = Bitset(size);
    nodeset.set();
    //initialise list of vertices to all vertices
    vertices = std::vector<int>(size);
    std::iota(vertices.begin(), vertices.end(), 0);
    //add edges to matrix
    for(int i = 0; i < m; i++){
        gmatrix[edge_list[2 * i]].set( edge_list[2 * i + 1]);
        gmatrix[edge_list[2 * i + 1]].set( edge_list[2 * i]);
    }

    bag.resize(size);
    vertex_rep.resize(size);
    for (int i = 0; i < size; ++i) {
        bag[i].push_back(i);
        vertex_rep[i] = i;
    }
    current_level = 0;
    removed_vertices = {};
    added_edges_trail = {};
    added_edges_limiter = {};
    vertex_rep_trail = {};
    bag_sizes_trail = {};
}


bool MGraph::has_edge(const int u, const int v) const {
    assert( u < size);
    assert( v < size);
    assert(gmatrix[u][v] == gmatrix[v][u]);
    return gmatrix[u][v];
}

void MGraph::add_edge(const int u, const int v) {
    assert(not has_edge(u,v));
    gmatrix[u].set(v);
    gmatrix[v].set(u);
    if (current_level > 0) {
        added_edges_trail.insert(added_edges_trail.end(), {u,v});
    }
}

void MGraph::remove_edge(const int u, const int v) {
    assert(has_edge(u,v));
    gmatrix[u].reset(v);
    gmatrix[v].reset(u);
}

bool MGraph::is_contracted(const int u, const int v) const {
    assert( u < size);
    assert( v < size);
    assert(not has_edge(u,v));
    return vertex_rep[u] == vertex_rep[v];
}

void MGraph::contract_vertices(int u, int v) {
    assert( not is_contracted(u, v));
    assert(vertex_rep[u] == u);
    assert(vertex_rep[v] == v);
    if (u > v){
        //want to contract into the smaller vertex
        std::swap(u,v);
    }

    //for w in N(v)\N(u), and u' in bag of u, add edge (u', w)
    nv_without_nu = setminus(gmatrix[v], gmatrix[u]);
    for (int uprime : bag[u]) {
        for (int w = nv_without_nu.find_first(); w != Bitset::npos; w = nv_without_nu.find_next(w)){
            add_edge(uprime, w);
        }
    }
    //for w in N(u)\N(v), and v' in bag of v, add edge (w, v')
    nu_without_nv = setminus(gmatrix[u], gmatrix[v]);
    for (int vprime : bag[v]) {
        for (int w = nu_without_nv.find_first(); w != Bitset::npos; w = nu_without_nv.find_next(w)){
            add_edge(vprime, w);
        }
    }

    assert(gmatrix[u] == (gmatrix[u] | gmatrix[v]));

    //update representing vertex and bags
    for (int vprime : bag[v]) {
        vertex_rep[vprime] = u;
    }
    std::copy(bag[v].begin(), bag[v].end(), std::back_inserter(bag[u]));

    nodeset.reset(v); //vertex contracted, not available anymore
    removed_vertices.push_back(v);
    num_vertices--;
    //contracted v into u, remove v from list of vertices
    auto eraseIt = std::remove(vertices.begin(), vertices.end(), v); //move vertex to end
    vertices.erase(eraseIt, vertices.end()); //erase vertex
}

void MGraph::separate_vertices(const int u, const int v) {
    if(has_edge(u,v)) {
        return;
    }
    //have to add edges also for other vertices in the bags
    for(int vp : bag[v]) {
        if(not has_edge(u, vp)) {
            add_edge(u, vp);
        }
    }
    for(int up : bag[u]) {
        if(not has_edge(v, up)) {
            add_edge(v, up);
        }
    }
}

double MGraph::density() const {
    int nvertices = vertices.size();
    double edges = 0.0;
    for(int vertex : vertices) {
        edges += (gmatrix[vertex]&nodeset).count(); //edges are already counted twice this way
    }
    return edges/(nvertices * (nvertices - 1));
}


//these backtracking eatures are inspired by the implementation of the original paper
void MGraph::notify_new_level() {
    current_level++;
    // std::cout << "notify new level " << current_level << "\n";
    assert(current_level == added_edges_limiter.size() + 1);
    assert(vertex_rep_trail.size() == bag_sizes_trail.size());

    //do not resize vertex_rep_trail etc. when backtracking so we might have already allocated all the space
    // otherwise, do allocate
    if(current_level > vertex_rep_trail.size()) {
        //assert it only became larger by one
        assert(current_level == vertex_rep_trail.size() + 1);
        vertex_rep_trail.emplace_back();
        bag_sizes_trail.emplace_back(size);
    }

    //new level means new limiter for added edges
    added_edges_limiter.push_back(static_cast<int>(added_edges_trail.size()));

    //copy current vertex representatives and bag sizes
    vertex_rep_trail[current_level - 1] = vertex_rep;
    for(auto w : vertices) {
        bag_sizes_trail[current_level - 1][w] = static_cast<int>(bag[w].size());
    }
    assert(vertex_rep_trail.size() == bag_sizes_trail.size());
    assert(added_edges_limiter.size() == current_level);
    assert(check_consistency());
}


void MGraph::notify_backtrack_level(const int new_level) {
    current_level = new_level;
    // std::cout << "Backtrack to level " << new_level << "\n";
    assert(new_level <= added_edges_limiter.size()); //can only backtrack on that level if the information is there

    // for each added edge until limiter, remove it again from the graph
    int edge_trail_size = added_edges_limiter[new_level];
    for (int i = edge_trail_size; i < added_edges_trail.size(); i += 2) {
        remove_edge(added_edges_trail[i], added_edges_trail[i+1]);
    }
    added_edges_trail.resize(edge_trail_size);
    added_edges_limiter.resize(new_level);

    //reset to copied versions of vertex_rep, and bags via bag_sizes_trail (no need to resize, see notify_new_level)
    vertex_rep = vertex_rep_trail[new_level];
    while (!removed_vertices.empty()) {
        int v = removed_vertices.back();
        assert(!nodeset[v]);
        if (vertex_rep[v] != v) {
            //only restore vertex if it was its own rep in the new_level
            break;
        }
        nodeset.set(v);
        removed_vertices.pop_back();
        num_vertices++;
        vertices.push_back(v); //restored vertex, add it back to list
    }
    //keep vector of vertices sorted
    std::sort(vertices.begin(), vertices.end());
    for(int w : vertices) {
        assert(vertex_rep[w] == w); //make sure active nodes are the representatives
        assert(bag_sizes_trail[new_level][w] > 0); //make sure the bag contains at least the vertex itself
        bag[w].resize( bag_sizes_trail[new_level][w]);
    }
    assert(check_consistency());
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


bool MGraph::check_consistency() {
    // std::cout << "checking graph consistency\n";
    for (int i = 0; i < size; ++i) {
        assert((!nodeset[i] or vertex_rep[i] == i)
            and (nodeset[i] or vertex_rep[i] != i));
    }

    int nvertices = nodeset.count();
    double edges = 0.0;
    for(int vertex = nodeset.find_first(); vertex != Bitset::npos; vertex = nodeset.find_next(vertex)) {
        edges += (gmatrix[vertex]&nodeset).count(); //edges are counted twice this way
    }
    assert(nvertices == num_vertices);

    Bitset bs(size);
    bs.set();
    for (int v : removed_vertices) {
        assert(!nodeset[v]);
        bs.reset(v);
    }
    assert(bs == nodeset);

    std::vector<int> nodeset_vertices;
    for(auto v = nodeset.find_first(); v!= Bitset::npos; v = nodeset.find_next(v)) {
        nodeset_vertices.push_back(v);
    }
    assert(nodeset_vertices == vertices);
    assert(nodeset.count() == vertices.size());

    for(auto v = nodeset.find_first(); v!= Bitset::npos; v = nodeset.find_next(v)) {
        assert(!bag[v].empty());
        assert(bag[v][0] == v);
        assert(vertex_rep[v] == v);
        bs.reset();
        for (int u : bag[v]) {
            assert(vertex_rep[u] == v);
            bs |= (gmatrix[u]);
        }
        if (!bs.is_subset_of(gmatrix[v])) {
            std::cout << "    bs[" << v << "] = " << bs << "\n";
            std::cout << "matrix[" << v << "] = " << gmatrix[v] << "\n";
        }
        assert(bs.is_subset_of(gmatrix[v]));
    }

    for(auto v = nodeset.find_first(); v!= Bitset::npos; v = nodeset.find_next(v)) {
        for(auto u = nodeset.find_first(); u!= Bitset::npos; u = nodeset.find_next(u)) {
            if (u == v) {
                continue;
            }
            if (gmatrix[v][u]) {
                assert(gmatrix[u][v]);
                for (int vp : bag[v]) {
                    if (!gmatrix[u][vp]) {
                        std::cout << "u = " << u << " v = " << v
                                  << " vp = " << vp << "\n";
                        std::cout << "matrix[v] = " << gmatrix[v] << '\n';
                        std::cout << "matrix[u] = " << gmatrix[u] << '\n';
                        std::cout << "bag[v] = " << bag[v] << '\n';
                        std::cout << "bag[u] = " << bag[u] << '\n';
                        assert(gmatrix[u][vp]);
                    }
                }
                for (int up : bag[u]) {
                    if (!gmatrix[v][up]) {
                        std::cout << "matrix[v] = " << gmatrix[v] << '\n';
                        std::cout << "matrix[u] = " << gmatrix[u] << '\n';
                        assert(gmatrix[v][up]);
                    }
                }
            }
        }
    }
    return true;
}

//compute setminus a\b = a intersect b^c
Bitset MGraph::setminus(const Bitset &a, const Bitset &b) {
    return a - (a & b);
}

void MGraph::print() {
    std::cout << "print matrix[i][j]\n";
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            std::cout << gmatrix[i][j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "Next, print vertex bags\n";
    for (const auto& vbag : bag) {
        for (int v : vbag) {
            std::cout << v << " ";
        }
        std::cout << "\n";
    }
}

int MGraph::dsatur_coloring(const Bitset &clique) const {
    std::vector<std::set<int>> coloring;
    Bitset available_vertices = nodeset & ~clique;
    std::vector<int> uncolored_vertices;
    for (int v = available_vertices.find_first(); v != Bitset::npos; v = available_vertices.find_next(v)) {
        uncolored_vertices.push_back(v);
    }
    int colors = uncolored_vertices.size() + clique.count();
    //initialise all vertices to be uncolored and have zero saturation
    std::vector<int> vertex_color(size, -1);
    std::vector<Bitset> forbidden_colors(colors, Bitset(size)); //colors rows and vertices columns!
    //record for an uncolored vertex to how many colored vertices it is connected to
    std::vector<int> saturation_level(size, 0);

    int nclique_colors = 0;
    if(clique.count() > 2){//ignore trivial or empty clique
        //fix color of the vertices of supplied clique and set saturation levels
        int color = 0;
        for(auto vertex = clique.find_first(); vertex != Bitset::npos; vertex = clique.find_next(vertex)){
            coloring.emplace_back(); //new color class
            coloring[color].insert(vertex); //assign first color to max degree vertex
            vertex_color[vertex] = color;
            forbidden_colors[color] |= gmatrix[vertex];
            color++;
            //update saturation levels
            for(auto neighbor = gmatrix[vertex].find_first(); neighbor != Bitset::npos; neighbor = gmatrix[vertex].find_next(neighbor)){
                saturation_level[neighbor]++;
            }
            nclique_colors++;
        }
    }
    //iterate while there are uncolored vertices
    while (not uncolored_vertices.empty()) {
        int max_saturated_vertex = 0;
        int max_saturated = -1;
        //find vertex with largest saturation degree
        for(int vertex : uncolored_vertices){
            //iterate over only uncolored vertices
            assert(vertex_color[vertex] == -1);
            if(saturation_level[vertex] > max_saturated){
                //found vertex with higher saturation, update
                max_saturated = saturation_level[vertex];
                max_saturated_vertex = vertex;
            } else if(saturation_level[vertex] == max_saturated and gmatrix[vertex].count() > gmatrix[max_saturated_vertex].count()){
                    //found vertex with equal saturation but higher degree, update
                    max_saturated_vertex = vertex;
            }
        }
        //find lowest color for vertex
        int lowest_color = -1;
        for(int color = 0; color < colors; color++) {
            if(forbidden_colors[color][max_saturated_vertex]) {
                continue; //color not available
            }
            lowest_color = color;
            break;
        }
        assert(lowest_color != -1);
        //assign color and update forbidden colors
        uncolored_vertices.erase(std::find(uncolored_vertices.begin(), uncolored_vertices.end(), max_saturated_vertex));
        available_vertices.reset(max_saturated_vertex);
        if(lowest_color >= coloring.size()) {
            assert(lowest_color == coloring.size());
            coloring.emplace_back(); //new color class
        }
        //assign first color to max saturated vertex
        coloring[lowest_color].insert(max_saturated_vertex);
        vertex_color[max_saturated_vertex] = lowest_color;
        forbidden_colors[lowest_color] |= gmatrix[max_saturated_vertex];

        //new vertex has been selected so update the saturation levels of all its neighbors if it's a different color
        //only iterate over neighbours that are still available
        Bitset neighbour_set = gmatrix[max_saturated_vertex] & available_vertices;
        for (auto vertex = neighbour_set.find_first(); vertex != Bitset::npos; vertex = neighbour_set.find_next(vertex)) {
            bool is_new_color = true;
            for (auto adj_vertex = gmatrix[vertex].find_first(); adj_vertex != Bitset::npos; adj_vertex = gmatrix[vertex].find_next(adj_vertex)) {
                if(adj_vertex != max_saturated_vertex and vertex_color[adj_vertex] == lowest_color){
                    //found neighbor of neighbor that uses msv color, so we don't increase the saturation level
                    is_new_color = false;
                    break;
                }
            }
            if(is_new_color) {
                saturation_level[vertex]++;
            }
        }
    }

    //colored all vertices
    return coloring.size();
}


int MGraph::sequential_coloring(const Bitset &clique) const {
    std::vector<std::set<int>> coloring;
    std::vector<int> available_vertices;
    available_vertices.reserve(nodeset.count());
    Bitset nodes_minus_clique = (nodeset & ~clique);
    assert(nodes_minus_clique.count() == (nodeset.count() - clique.count()));
    for(int vertex = nodes_minus_clique.find_first(); vertex != Bitset::npos; vertex = nodes_minus_clique.find_next(vertex)) {
        available_vertices.push_back(vertex);
    }
    //sort by decreasing bag size and then by degree
    std::sort(available_vertices.begin(), available_vertices.end(),
        [&](int v, int w){return bag[v].size() > bag[w].size() or (bag[v].size() == bag[w].size() and gmatrix[v].count() > gmatrix[w].count());});

    int colors = available_vertices.size() + clique.count();
    //initialise all vertices to be uncolored and have zero saturation
    std::vector<int> vertex_color(size, -1);
    std::vector<Bitset> forbidden_colors(colors, Bitset(size)); //colors rows and vertices columns!

    int nclique_colors = 0;
    if(clique.count() > 2){//ignore trivial or empty clique
        //fix color of the vertices of supplied clique and set saturation levels
        int color = 0;
        for(auto vertex = clique.find_first(); vertex != Bitset::npos; vertex = clique.find_next(vertex)){
            coloring.emplace_back(); //new color class
            coloring[color].insert(vertex); //assign first color to max degree vertex
            vertex_color[vertex] = color;
            forbidden_colors[color] |= gmatrix[vertex];
            color++;
            nclique_colors++;
        }
    }

    //iterate over all uncolored vertices
    for(int vertex : available_vertices){
        assert(vertex_color[vertex] == -1);
        //find lowest color for vertex
        int lowest_color = -1;
        for(int color = 0; color < colors; color++) {
            if(forbidden_colors[color][vertex]) {
                continue; //color not available
            }
            lowest_color = color;
            break;
        }
        assert(lowest_color != -1);
        //assign color and update forbidden colors
        if(lowest_color >= coloring.size()) {
            assert(lowest_color == coloring.size());
            coloring.emplace_back(); //new color class
        }
        coloring[lowest_color].insert(vertex); //assign first color to max degree vertex
        vertex_color[vertex] = lowest_color;
        forbidden_colors[lowest_color] |= gmatrix[vertex];
    }
    //colored all vertices
    return coloring.size();
}

int MGraph::IS_extract(const Bitset &clique) const {
    //try to build a good coloring by building large independent sets
    //iterate over vertices and try to add it to first independent set that admits it, or create a new one if we have to
    //kind of the complement of the greedy clique search
    std::vector<Bitset> independent_sets;

    std::vector<int> order = vertices; //generate lexikographic order on remaining vertices, later other heuristic
    assert(not order.empty());
    //sort in order of decreasing bag size
    // std::sort(order.begin(), order.end(), [&](int v, int w){return bag[v].size() > bag[w].size();});
    // std::sort(order.begin(), order.end(), [&](int v, int w){return gmatrix[v].count() > gmatrix[w].count();});

    //create singleton independent sets for the clique vertices
    if(clique.count() > 2){//ignore trivial or empty clique
        for(auto vertex = clique.find_first(); vertex != Bitset::npos; vertex = clique.find_next(vertex)){
            independent_sets.emplace_back(size);
            independent_sets.back().set(vertex);
        }
    }

    //iterate over vertices and add to independent set if possible
    for (int v : order) {
        if(clique[v]) {
            continue;
        }
        bool inserted = false;
        for(auto & indset : independent_sets){
            //check whether vertex is adjacent to none of the vertices, i.e. N(v) intersect indset == empty
            if((gmatrix[v] & indset).none()) {
                indset.set(v);
                inserted = true;
                break;
            }
        }
        //if not possible to add vertex to existing independent set, create new one or check if we can recolor
        if(not inserted) {
            bool recolored = try_recolor(v, independent_sets);
            if(not recolored) {
                independent_sets.emplace_back(size);
                independent_sets.back().set(v);
            }
        }
    }
    //no code to insert vertices into previous sets necessary here
    return independent_sets.size();
}

int MGraph::ISEQ() const {
    int color_index = 0;
    std::vector<Bitset> coloring; //store colorclasses as bitsets for now, convert them later
    Bitset nodes = nodeset; //copy available vertices
    Bitset tmp_nodes = nodes; //util set
    while (nodes.any()) { //while there are still uncolored vertices
        coloring.emplace_back(size);
        assert(color_index + 1 == coloring.size());
        while (tmp_nodes.any()) {
            int v = tmp_nodes.find_first();
            bool recolored = false;
            if(tmp_nodes.find_next(v) == Bitset::npos) {
                //v is only vertex left
                recolored = try_recolor(v, coloring);
            }
            tmp_nodes.reset(v);
            if(not recolored) {
                coloring[color_index].set(v); //assign color
                tmp_nodes &= ~gmatrix[v]; //remove adjacent vertices
            }
            nodes.reset(v);
        }
        //made colorclass as large as possible, move on to next with available vertices
        tmp_nodes = nodes;
        color_index++;
    }
    //nodes is empty, every vertex is colored and coloring.size() colors were used
    return coloring.size();
}

bool MGraph::try_recolor(const int vertex, std::vector<Bitset> &coloring) const {
    //last class of coloring is the new colorclass for vertex
    for (auto k1 = coloring.begin(); k1 != coloring.end() - 2; ++k1) {
        auto ck_intersect_nv = *k1 & gmatrix[vertex];
        int w = ck_intersect_nv.find_first();
        //check that w is Only element of intersection
        if(w != Bitset::npos and ck_intersect_nv.find_next(w) == Bitset::npos) {
            for (auto k2 = std::next(k1); k2 != coloring.end() - 1; ++k2) {
                if((*k2 & gmatrix[w]).none()) {
                    k1->reset(w);
                    k1->set(vertex);
                    k2->set(w);
                    return true;
                    //called double swap
                }
            }
        }
        //intersection was empty, does normally not happen but can be caused by previous swap above in an earlier call
        else if(w == Bitset::npos) {
            k1->set(vertex);
            return true;
            //called single swap
        }
    }
    return false;
}

int MGraph::greedy_cliques(std::vector<Bitset> &clique_list, const int clique_limit) {
    //we maintain a list of cliques and iterate through vertices in given ordering
    //try to add vertex to all cliques in list or create new one if it can't be added

    std::vector<Bitset> cliques;
    std::vector<int> clique_sz; //to keep track of sizes instead of needing to call Bitset.count() which can be linear
    std::vector<int> last_inserted_index(size); //track the last clique where this vertex was added

    assert(not vertices.empty());
    // std::vector<int> order = vertices;
    // sorting not used currently as this would cause quite the overhead since cliques are computed so often
    //sort in order of decreasing bag size or degree
    // std::sort(order.begin(), order.end(), [&](int v, int w){return mgraph.bag[v].size() > mgraph.bag[w].size();});
    // std::sort(order.begin(), order.end(), [&](int v, int w){return mgraph.gmatrix[v].count() > mgraph.gmatrix[w].count();});

    //iterate over vertices and add to cliques if possible
    for (int v : vertices) {
        bool inserted = false;
        for(int i = 0; i < cliques.size(); i++){
            Bitset &clq = cliques[i];
            //check whether vertex is adjacent to all vertices of clique, i.e. N(v) intersect clq == clq
            if((gmatrix[v] & clq) == clq) {
                clq.set(v);
                clique_sz[i]++;
                last_inserted_index[v] = i;
                inserted = true;
            }
        }
        //if not possible to add vertex to existing clique, create new one
        if(not inserted and cliques.size() < clique_limit) {
            last_inserted_index[v] = static_cast<int>(cliques.size());
            cliques.emplace_back(size);
            cliques.back().set(v);
            clique_sz.push_back(1);
        }
    }
    //iterate over vertices in order once more, since vertices were not attempted to be added to cliques created later
    for (int v : vertices) {
        for (int i = last_inserted_index[v]; i < cliques.size(); i++) {
            Bitset &clq = cliques[i];
            if((gmatrix[v] & clq) == clq) {
                clq.set(v);
                clique_sz[i]++;
                last_inserted_index[v] = i;
            }
        }
    }

    //store all largest cliques that were found in clique_list
    int max_size = std::distance(clique_sz.begin(), std::max_element(clique_sz.begin(), clique_sz.end()));
    assert(cliques[max_size].count() == clique_sz[max_size]);
    for (int i = max_size; i < cliques.size(); i++) {
        assert(is_clique(cliques[i]));
        if(clique_sz[i] == clique_sz[max_size]) {
            clique_list.push_back(cliques[i]);
        }
    }
    return clique_sz[max_size];
}

bool MGraph::is_clique(const Bitset &clique) const {
    assert(clique.size() == size);
    for (int v = clique.find_first(); v != Bitset::npos; v = clique.find_next(v)) {
        for (int w = clique.find_next(v); w != Bitset::npos; w = clique.find_next(w)) {
            if (not has_edge(v,w)) {
                return false;
            }
        }
    }
    return true;
}

MGraph::SubGraph::SubGraph(const Bitset &clique)
    : SubGraph(clique.size()){
    for(auto vertex = clique.find_first(); vertex != Bitset::npos; vertex = clique.find_next(vertex)) {
        matrix[vertex] |= clique; //add all neighbours, i.e. all clique vertices
        matrix[vertex].reset(vertex); //remove vertex, not neighbour of itself
        num_vertices++; //added a node
        nodes.push_back(vertex);
        nodeset.set(vertex);
    }
}


int MGraph::mycielsky_extension_clique(SubGraph &subgraph, const int threshold) {
    //function to extend an initial clique to larger generalised mycielsky graph
    // returned is the number of successful iterations and stores reason clause in external_reasons if succesful
    int iterations = 0;

    Bitset W = Bitset(size);
    while (subgraph.num_vertices < size and iterations < threshold){ //can stop when enough iterations were succesful
        //initialize W=V and S_v = {v}
        W.set();
        std::vector<Bitset> vecS(subgraph.num_vertices, Bitset(size));
        for(int i = 0; i < subgraph.num_vertices; i++){
            int v = subgraph.nodes[i];
            vecS[i].set(v);
        }
        //first loop computing S_v and intersecting W
        for(int i = 0; i < subgraph.num_vertices; i++){
            int v = subgraph.nodes[i];
            Bitset NG_Sv = gmatrix[v]; //so we can skip v in next loop
            for (int u : vertices) {
                if(u == v) {
                    continue;
                }
                if (subgraph.matrix[v].is_subset_of(gmatrix[u])) {
                    vecS[i].set(u);
                    NG_Sv |= gmatrix[u]; //union
                }
            }
            W &= NG_Sv; //intersection
            if(W.none()){ //W already empty, can stop here
                break;
            }
        }
        //case that W is not empty, can extend subgraph and increase lower bound
        if(W.any()){
            iterations++;
            //instead of copying the subgraph, we work with current one and remember which vertices/edges to add and do so later
            std::vector<int> new_vertices;
            std::vector<std::pair<int, int>> new_edges;
            int w = W.find_first();
            if(not subgraph.nodeset[w]) {
                new_vertices.push_back(w);
            }

            for(int i = 0; i < subgraph.num_vertices; i++){
                int v = subgraph.nodes[i];
                //insert vertex u from the intersection
                Bitset  NG_w_intersect_S_v = gmatrix[w] & vecS[i];
                int u = NG_w_intersect_S_v.find_first();
                if(not subgraph.nodeset[u]) {
                    new_vertices.push_back(u);
                }
                //add edges between u and w, and u and N_H(v)
                new_edges.push_back({u,w});
                for(int vn = subgraph.matrix[v].find_first(); vn != Bitset::npos; vn = subgraph.matrix[v].find_next(vn)){
                    new_edges.emplace_back(u,vn);
                }
            }

            //collectecd vertices and edges to add, now add them to subgraph and finally sort them
            for (int u : new_vertices) {
                subgraph.nodes.push_back(u);
                subgraph.num_vertices++;
                subgraph.nodeset.set(u);
            }
            for (auto [u, vn] : new_edges) {
                subgraph.matrix[u].set(vn);
                subgraph.matrix[vn].set(u);
            }
            std::sort(subgraph.nodes.begin(), subgraph.nodes.end());
        }
        else {
            break;
        }
    }
    assert(iterations <= threshold);
    return iterations;
}

std::vector<Bitset> MGraph::get_active_subgraph(std::vector<int> &index_mapping) const {
    index_mapping.clear();
    index_mapping.resize(num_vertices);
    std::vector<Bitset> active_subgraph(num_vertices, Bitset(num_vertices)); //n by n matrix on active vertices
    for(int i = 0; i < num_vertices; i++) {
        int active_v = vertices[i];
        index_mapping[i] = active_v;
        for(int j = i + 1; j < num_vertices; j++) {
            int active_u = vertices[j];
            if (gmatrix[active_v][active_u]) {
                active_subgraph[i].set(j);
                active_subgraph[j].set(i);
            }
        }
    }
    return active_subgraph;
}

int MGraph::mnts_clique(std::vector<Bitset> &clique_list, int aim, int mnts_length, int mnts_depth, int random_seed) const {
    std::vector<int> index_mapping;
    std::vector<Bitset> active_subgraph = get_active_subgraph(index_mapping);
    MNTS mnts(active_subgraph, aim, mnts_length, mnts_depth, random_seed);
    Bitset clique = mnts.runSearch();
    // map clique from subgraph indices back to the original graph
    Bitset clique_original_graph(size);
    for (int u = clique.find_first(); u != Bitset::npos; u = clique.find_next(u)) {
        clique_original_graph.set(index_mapping[u]);
    }
    assert(is_clique(clique_original_graph));
    clique_list = {clique_original_graph};
    return clique_original_graph.count();
}




