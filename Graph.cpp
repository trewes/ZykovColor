#include "Graph.h"


namespace Graph{

template<typename T>
std::ostream &operator<<(std::ostream &s, const std::set<T> &set) {
    s << "{";
    std::string sep;
    for(T el : set){
        s << sep << el;
        sep = ", ";
    }
    return s << "}";
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
T random_element(const std::vector<T> &v) {
    auto it = v.cbegin();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(0, v.size() - 1);
    int random = distr(gen);
    std::advance(it, random);
    return *it;
}

Permutation identity(int size) {
    Permutation identity(size);
    for(int i = 0; i < size; i++)
        identity[i] = i;
    return identity;
}

Permutation perm_inverse(const Permutation &perm) {
    Permutation inverse(perm.size());
    for(int i = 0; i < static_cast<int>(perm.size()); i++){
        inverse[perm[i]] = i;
    }
    return inverse;
}

std::vector<VertexType> apply_permutation(std::vector<VertexType> vertices, Permutation perm) {
    for(VertexType & vertex : vertices){
        vertex = perm[vertex];
    }
    return vertices;
}


NeighborList Graph::get_neighbor_list() const {
    NeighborList list;
    list.resize(_ncount);
    for(VertexType i = 0; i < _ecount; i++){
        list[_elist[2 * i]].push_back( _elist[2 * i + 1]);
        list[_elist[2 * i + 1]].push_back( _elist[2 * i]);
    }
    for(VertexType i = 0; i < _ncount; i++){//like it better if the neighbors are sorted
        std::sort(list[i].begin(), list[i].end());
    }
    return list;
}

Graph::Graph(VertexType vertex_count, VertexType edge_count, std::vector<VertexType> edge_list) :
        _ncount(vertex_count), _ecount(edge_count), _elist(std::move(edge_list)) {}

VertexType Graph::ncount() const {
    return _ncount;
}

VertexType Graph::ecount() const {
    return _ecount;
}

double Graph::density() const {
    return static_cast<double>(2*_ecount)/(_ncount *(_ncount - 1));
}

void Graph::print() const {
    std::cout << "Graph has " << _ncount << " vertices and " << _ecount << " edges.\n";
    NeighborList neighbors = get_neighbor_list();
    for(auto &adj : neighbors){
        std::cout << adj << "\n";
    }
}

std::vector<VertexType> Graph::elist() const {
    return _elist;
}

std::vector<VertexType> Graph::complement_elist() const {
    NeighborList neighbors = get_neighbor_list();
    VertexType complement_ecount = (_ncount * (_ncount - 1)) / 2 - _ecount;
    std::vector<VertexType> complement_edge_list(2 * complement_ecount);
    VertexType cecount = 0;//count edges we add
    //iterate over all pairs of vertices
    for(VertexType v = 0; v < _ncount; v++){
        VertexType a = 0; //neighbor vertex
        VertexType a_i  = 0; //index of neighbor
        for (VertexType w = v + 1; w < _ncount; ++ w) {
            //find first neighbor that is greater equal to w
            while (a_i < static_cast<int>(neighbors[v].size()) && a < w) {
                a = neighbors[v][a_i];
                ++a_i;
            }
            //if neighbor 'a' is equal to w, w is adjacent to v
            //if not, insert edge in complementary graph
            if (w != a) {
                complement_edge_list[2*cecount]   = v;
                complement_edge_list[2*cecount+1] = w;
                ++cecount;
            }
        }
    }
    if(cecount != complement_ecount){
        throw std::runtime_error("Did not produce the right amount of complement edges.");
    }
    return complement_edge_list;
}

Graph Graph::get_complement() const {
    //get all edges not in this graph and build the complement edge list from that
    //use that to return the complement graph
    VertexType complement_ecount = (_ncount * (_ncount - 1)) / 2 - _ecount;
    return {_ncount, complement_ecount,  complement_elist()};
}


Graph::Graph(const char *filename) : _ncount(0), _ecount(0) {
    std::string string_name = filename;
    if(string_name.substr(string_name.find_last_of('.') + 1) == "col"){
        read_dimacs(filename);
    } else if(string_name.substr(string_name.find_last_of('.') + 1) == "g6"){
        read_graph6(filename);
    } else if(string_name.substr(string_name.find_last_of('.') + 1) == "b"){
        throw std::runtime_error("Binary dimacs format is currently not supported, please translate it yourself.");
    } else{
        std::cout << "c Warning : Filename does not have a recognized ending, trying to read as DIMACS graph!\n";
        read_dimacs(filename);
    }
    simplify_graph();
}

Graph::Graph(std::string g6_string) : _ncount(0), _ecount(0) {
    read_graph6_string(std::move(g6_string));
}


void Graph::read_dimacs(const char *filename) {
    std::ifstream file(filename);
    if(not file){
        throw std::runtime_error("Cannot open file.");
    }

    std::string line;
    std::getline(file, line);//read in the file line by line

    while(line[0] == 'c' or line[0] == '\r' or line[0] == '\0'){          //special character to account for empty lines
        std::getline(file, line);                                                               //ignore comments
    }
    if(line[0] != 'p'){
        throw std::runtime_error("File is not in correct dimacs format.");          //a different case was expected
    }
    std::string p, edge;
    int n, e;
    std::stringstream ss(line);
    ss >> p >> edge >> n >> e;                                                       //read in number of nodes and edges

    _ncount = n;
    _ecount = e;
    _elist.reserve(2 * _ecount);

    std::getline(file, line);
    while(line[0] == 'n'){
        std::cout << "The dimacs file is specifying colors of vertices, these will be ignored for this problem.\n";
        std::getline(file, line);
    }

    char first;
    VertexType head, tail;
    do{
        while(line[0] == 'c' or line[0] == '\r' or
              line[0] == '\0'){          //special character to account for empty lines
            std::getline(file, line);                                                           //ignore comments
        }
        ss.clear();
        ss.str(std::string());
        ss << line;
        ss >> first >> head >> tail;
        _elist.push_back(head - 1);
        _elist.push_back(tail - 1);
    } while(std::getline(file, line));
}

void Graph::read_graph6(const char *filename) {
    std::ifstream file(filename);
    if(not file){
        throw std::runtime_error("Cannot open file.");
    }
    std::cout << "Warning: if the graph6 file contains more than one graph, the first is read in\n";
    std::string line;
    std::getline(file, line);
    read_graph6_string(line);
}

//this function of reading in the graph6 format was taken from treedecomposition.com
void Graph::read_graph6_string(std::string g6_string) {// Extract vertex count from graph6
    if(g6_string.substr(0, 1) == ":"){
        g6_string = g6_string.substr(1, g6_string.size());
    }

    int vertexCount = -1;

    int r0 = g6_string[0];
    int adjIdx = 1;

    // single byte expansion (0 <= n <= 62)
    if(r0 >= 0 + 63 and r0 <= 63 + 62){
        vertexCount = r0 - 63;
        adjIdx = 1;
    }
        // four or eight byte expansions
    else if(r0 == 126){
        char r1 = g6_string[1];

        // eight byte expansions (258048 <= n <= 68719476735)
        if(r1 == 126){
            adjIdx = 8;
            throw std::runtime_error("Terribly sorry, but we don\t do 8-byte expansions yet :(");
        }
            // four byte expansion (63 <= n <= 258047)
        else{
            if(g6_string.size() < 4){
                throw std::runtime_error("This input seems to be too short for a 4-byte expansion :(");
            } else{
                vertexCount = 0;
                adjIdx = 4;

                for(int i = 0; i < 3; i++){
                    vertexCount = vertexCount | ((g6_string[1 + i] - 63) << ((2 - i) * 6));
                }
            }
        }
    }

    // Only start working when we actually have a vertex count
    if(vertexCount > -1){
        _ncount = vertexCount;
        int edgeCount = 0;
        // Extract adjacency list from graph6,
        for(int n = 1; n < vertexCount; n++){
            // Natural numbers sum n*(n+1)/2, but only need previous n so (n-1)*(n+1-1)/2
            int nPos = n * (n - 1) / 2;
            for(int m = 0; m < n; m++){
                int bitPos = nPos + m;
                int charIndex = adjIdx + int(std::floor(bitPos / 6));
                int charVal = g6_string[charIndex] - 63;

                int bitIndex = bitPos % 6;
                int bitShift = (5 - bitIndex);
                int bitMask = (1 << bitShift);

                int bitVal = ((charVal & bitMask) >> bitShift);

                // If the bit for this particular position is 1 then {m, n} is an edge
                if(bitVal == 1){
                    _elist.push_back(n);
                    _elist.push_back(m);
                    edgeCount++;
                }
            }
        }
        _ecount = edgeCount;
    } else{
        throw std::runtime_error("Was not able to get vertex count from graph6 file.");
    }
}


void Graph::simplify_graph(){

    NeighborList neighbors = get_neighbor_list();
    std::vector<VertexType> simplified_elist;
    simplified_elist.reserve(2*_ecount);
    int  simplified_ecount = 0;
    for(VertexType v = 0; v < _ncount; v++){
        if(neighbors[v].empty())
            continue;
        //add first edge if it is not a loop and edge from smaller to larger vertex
        if(v < neighbors[v][0]){
            simplified_elist.push_back(v);
            simplified_elist.push_back(neighbors[v][0]);
            simplified_ecount++;
        }
        for(VertexType j = 1; j < static_cast<int>(neighbors[v].size()); j++ ){
            //test for loops or duplicate edges, uses that the adjacency lists from get_neighbor_list are sorted
            //only add edge from smalle to larger vertex, also skip otherwise
            if( (v < neighbors[v][j] ) and (neighbors[v][j - 1] != neighbors[v][j]) ){
                simplified_elist.push_back(v);
                simplified_elist.push_back(neighbors[v][j]);
                simplified_ecount++;
            }
        }
    }
    if(simplified_ecount != (static_cast<int>(simplified_elist.size()) / 2)){
        throw std::runtime_error("Did not simplify the graph correctly.");
    }
    _elist = simplified_elist;
    _elist.shrink_to_fit();
    _ecount = simplified_ecount;
}


Graph Graph::perm_graph(const Permutation &perm) const {
    if(static_cast<VertexType>(perm.size()) != _ncount){
        throw std::runtime_error("Size of graph and permutation do not match.");
    }
    Graph permuted_graph(_ncount, _ecount, {});
    permuted_graph._elist.resize(2 * _ecount);
    for(VertexType i = 0; i < 2 * _ecount; i++){
        permuted_graph._elist[i] = perm[_elist[i]];
    }

    //also update information about removed vertices/edges and vertex mapping
    //does not reflect the action being performed on removed vertices/edges
    //so before calling recover_reductions, we have to revert the permutation and get the original graph!
    permuted_graph.removed_vertices = removed_vertices;
    permuted_graph.removed_edges = removed_edges;
    permuted_graph.remap = remap;

    return permuted_graph;
}

Coloring Graph::dsatur(Permutation &ordering, const std::vector<VertexType> &clique) const {
    NeighborList neighbors = get_neighbor_list();
    return dsatur(ordering, neighbors, clique);
}

Coloring Graph::dsatur(Permutation &ordering, const NeighborList &neighbors, const std::vector<VertexType> &clique) const {
    ordering.clear();
    ordering.reserve(_ncount);
    Coloring coloring;
    //initialise all vertices to be uncolored and have zero saturation
    std::vector<int> vertex_color(_ncount, -1);
    //record for an uncolored vertex to how many colored vertices it is connected to
    std::vector<int> saturation_level(_ncount, 0);

    int nclique_colors = 0;
    if(clique.size() > 2){//ignore trivial or empty clique
        //fix color of the vertices of supplied clique and set saturation levels
        int color_index = 0;
        for(VertexType v : clique){
            ordering.push_back(v);
            coloring.emplace_back(); //new color class
            coloring[color_index].insert(v); //assign first color to max degree vertex
            vertex_color[v] = color_index;
            color_index++;
            //update saturation levels
            for(VertexType neighbor : neighbors[v]){
                saturation_level[neighbor]++;
            }
            saturation_level[v] = std::numeric_limits<int>::min();
        }
        nclique_colors = int(clique.size());
    }

    //main loop looking for the next vertex according to its saturation level
   for(int iteration = 0; iteration < _ncount - nclique_colors; iteration++){
        VertexType max_saturated_vertex = 0;
        int max_saturated = -1;
        //keep set of tied vertices to later choose one at random
        std::vector<VertexType> saturated_tied_max_degrees;
        for(VertexType vertex = 0; vertex < _ncount; vertex++){
            if(vertex_color[vertex] != -1)
                continue;

            if(saturation_level[vertex] > max_saturated){
                //found vertex with higher saturation, update
                max_saturated = saturation_level[vertex];
                max_saturated_vertex = vertex;

                saturated_tied_max_degrees.clear();
                saturated_tied_max_degrees.push_back(vertex);
            } else if(saturation_level[vertex] == max_saturated){
                if(neighbors[vertex].size() > neighbors[max_saturated_vertex].size()){
                    //found vertex with equal saturation but higher degree, update
                    max_saturated_vertex = vertex;

                    saturated_tied_max_degrees.clear();
                    saturated_tied_max_degrees.push_back(vertex);
                } else if(neighbors[vertex].size() == neighbors[max_saturated_vertex].size()){
                    //same saturation level, same degree, add to set of possible candidates
                    saturated_tied_max_degrees.push_back(vertex);
                }
            }
        }
        //select max saturated vertex, at random if there exist multiple
        if(saturated_tied_max_degrees.size() > 1 and random_tiebreaks){
            max_saturated_vertex = random_element(saturated_tied_max_degrees);
        }
        ordering.push_back(max_saturated_vertex);

       int lowest_color = 0;
       bool failed;
       do {
           failed = false;
           for (VertexType neighbor : neighbors[max_saturated_vertex]) {
               if (vertex_color[neighbor] == lowest_color) {
                   ++lowest_color;
                   failed = true;
                   break;
               }
           }
       } while (failed);


        if(lowest_color < int(coloring.size())){
            //we can add chosen vertex to some existing color class
            coloring[lowest_color].insert(max_saturated_vertex);
            vertex_color[max_saturated_vertex] = lowest_color;
        } else if(lowest_color == int(coloring.size())){
            //new color class is created for max_saturated_vertex
            //before that we try to recolor two vertices to not have to use a new color
            bool recolored = try_color_swap(max_saturated_vertex, neighbors, coloring, vertex_color);
            if(not recolored){
                //not able to recolor, add new color class
                coloring.emplace_back();
                coloring[lowest_color].insert(max_saturated_vertex);
                vertex_color[max_saturated_vertex] = lowest_color;
            }
        } else{
            std::cout << "Error: lowest color is too high, something must have gone wrong.\n";
            throw std::runtime_error("Error in Dsatur coloring");
        }

        //new vertex has been selected so update the saturation levels of all its neighbors if it's a different color
        int msv_color = vertex_color[max_saturated_vertex];
        for(VertexType vertex : neighbors[max_saturated_vertex]){
            bool is_new_color = true;
            //possibly check that neighbor has been selected, i.e. it's level is equal to numeric_limits<int>::min()
            if(saturation_level[vertex] != std::numeric_limits<int>::min()){
                for(VertexType adj_vertex : neighbors[vertex]){
                    if(adj_vertex != max_saturated_vertex and vertex_color[adj_vertex] == msv_color){
                        //found neighbor of neighbor that uses msv color, so we don't increase the saturation level
                        is_new_color = false;
                        break;
                    }
                }
                if(is_new_color)
                    saturation_level[vertex]++;
            }
        }
        //remove selected vertex from any further consideration
        saturation_level[max_saturated_vertex] = std::numeric_limits<int>::min();
    }
    //colored all vertices
    return coloring;
}


Coloring Graph::dsatur_original(Permutation &ordering, const std::vector<VertexType> &clique) const {
    NeighborList neighbors = get_neighbor_list();
    return dsatur_original(ordering, neighbors, clique);
}


Coloring
Graph::dsatur_original(Permutation &ordering, const NeighborList &neighbors, const std::vector<VertexType> &clique) const {
    ordering.clear();
    ordering.reserve(_ncount);
    Coloring coloring;
    //initialise all vertices to be uncolored
    std::vector<int> vertex_color(_ncount, -1);
    //record for an uncolored vertex to how many colored vertices it is connected
    //also track degree of degree in uncolored subgraph
    std::vector<int> saturation_level(_ncount, 0);
    std::vector<int> uncolored_subgraph_degree(_ncount, 0);
    for(VertexType v = 0; v < _ncount; v++)
        uncolored_subgraph_degree[v] = int(neighbors[v].size());

    int nclique_colors = 0;
    if(clique.size() > 2){//ignore trivial or empty clique
        //fix color of the vertices of supplied clique and set saturation levels
        int color_index = 0;
        for(VertexType v : clique){
            ordering.push_back(v);
            coloring.emplace_back(); //new color class
            coloring[color_index].insert(v); //assign first color to max degree vertex
            vertex_color[v] = color_index;
            color_index++;
            //update saturation levels
            for(VertexType neighbor : neighbors[v]){
                saturation_level[neighbor]++;
                uncolored_subgraph_degree[neighbor]--;
            }
            saturation_level[v] = std::numeric_limits<int>::min();
            uncolored_subgraph_degree[v] = std::numeric_limits<int>::min();
        }
        nclique_colors = int(clique.size());
    }


    //main loop looking for the next vertex according to its saturation level
    for(int iteration = 0; iteration < _ncount - nclique_colors; iteration++){
        VertexType max_saturated_vertex = 0;
        int max_saturated = -1;
        //keep set of tied vertices to later choose one at random
        std::vector<VertexType> saturated_tied_max_degrees;
        for(VertexType vertex = 0; vertex < _ncount; vertex++){
            if(vertex_color[vertex] != -1)
                continue;
            if(saturation_level[vertex] > max_saturated){
                //found vertex with higher saturation, update
                max_saturated = saturation_level[vertex];
                max_saturated_vertex = vertex;

                saturated_tied_max_degrees.clear();
                saturated_tied_max_degrees.push_back(vertex);
            } else if(saturation_level[vertex] == max_saturated){
                if(uncolored_subgraph_degree[vertex] > uncolored_subgraph_degree[max_saturated_vertex]){
                    //found vertex with equal saturation but higher degree, update
                    max_saturated_vertex = vertex;

                    saturated_tied_max_degrees.clear();
                    saturated_tied_max_degrees.push_back(vertex);
                } else if(neighbors[vertex].size() == neighbors[max_saturated_vertex].size()){
                    //same saturation level, same degree, add to set of possible candidates
                    saturated_tied_max_degrees.push_back(vertex);
                }
            }
        }
        //select max saturated vertex, at random if there exist multiple
        if(saturated_tied_max_degrees.size() > 1 and random_tiebreaks){
            max_saturated_vertex = random_element(saturated_tied_max_degrees);
        }
        ordering.push_back(max_saturated_vertex);

        int lowest_color = 0;
        bool failed;
        do {
            failed = false;
            for (VertexType neighbor : neighbors[max_saturated_vertex]) {
                if (vertex_color[neighbor] == lowest_color) {
                    ++lowest_color;
                    failed = true;
                    break;
                }
            }
        } while (failed);

        if(lowest_color < int(coloring.size())){
            //we can add chosen vertex to some existing color class
            coloring[lowest_color].insert(max_saturated_vertex);
            vertex_color[max_saturated_vertex] = lowest_color;
        } else if(lowest_color == int(coloring.size())){
            //new color class is created for max_saturated_vertex
            //before that we try to recolor two vertices to not have to use a new color
            bool recolored = try_color_swap(max_saturated_vertex, neighbors, coloring, vertex_color);
            if(not recolored){
                //not able to recolor, add new color class
                coloring.emplace_back();
                coloring[lowest_color].insert(max_saturated_vertex);
                vertex_color[max_saturated_vertex] = lowest_color;
            }
        } else{
            std::cout << "Error: lowest color is too high, something must have gone wrong.\n";
            throw std::runtime_error("Error in Dsatur coloring");
        }

        //new vertex has been selected so update the saturation levels of all its neighbors if it's a different color
        int msv_color = vertex_color[max_saturated_vertex];
        for(VertexType vertex : neighbors[max_saturated_vertex]){
            bool is_new_color = true;
            //possibly check that neighbor has been selected, i.e. it's level is equal to numeric_limits<int>::min()
            if(saturation_level[vertex] != std::numeric_limits<int>::min()){
                for(VertexType adj_vertex : neighbors[vertex]){
                    if(adj_vertex != max_saturated_vertex and vertex_color[adj_vertex] == msv_color){
                        //found neighbor of neighbor that uses msv color, so we don't increase the saturation level
                        is_new_color = false;
                        break;
                    }
                }
                if(is_new_color){
                    saturation_level[vertex]++;
                }
                //also update degree of neighbor in uncolored subgraph
                uncolored_subgraph_degree[vertex]--;
            }
        }
        //remove selected vertex from any further consideration
        saturation_level[max_saturated_vertex] = std::numeric_limits<int>::min();
        uncolored_subgraph_degree[max_saturated_vertex] = std::numeric_limits<int>::min();
    }
    //colored all vertices
    return coloring;
}


Coloring Graph::max_connected_degree_coloring(Permutation &ordering, const std::vector<VertexType> &clique) const {
    NeighborList neighbors = get_neighbor_list();
    return max_connected_degree_coloring(ordering, neighbors, clique);
}

Coloring Graph::max_connected_degree_coloring(Permutation &ordering, const NeighborList &neighbors,
                                              const std::vector<VertexType> &clique) const {
    ordering.clear();
    ordering.reserve(_ncount);
    Coloring coloring;
    //initialise all vertices to be uncolored
    std::vector<int> vertex_color(_ncount, -1);
    //record for an uncolored vertex to how many colored vertices it is connected
    std::vector<int> num_selected_neighbors(_ncount, 0);

    int nclique_colors = 0;
    if(clique.size() > 2){//ignore trivial or empty clique
        //fix color of the vertices of supplied clique and set saturation levels
        int color_index = 0;
        for(VertexType v : clique){
            ordering.push_back(v);
            coloring.emplace_back(); //new color class
            coloring[color_index].insert(v); //assign first color to max degree vertex
            vertex_color[v] = color_index;
            color_index++;
            //update saturation levels
            for(VertexType neighbor : neighbors[v]){
                num_selected_neighbors[neighbor]++;
            }
            num_selected_neighbors[v] = std::numeric_limits<int>::min();
        }
        nclique_colors = int(clique.size());
    }


    //main loop looking for the next vertex according to its saturation level
    for(int iteration = 0; iteration < _ncount - nclique_colors; iteration++){
        VertexType max_saturated_vertex = 0;
        int max_saturated = -1;
        //keep set of tied vertices to later choose one at random
        std::vector<VertexType> saturated_tied_max_degrees;
        for(VertexType vertex = 0; vertex < _ncount; vertex++){
            if(vertex_color[vertex] != -1)
                continue;

            if(num_selected_neighbors[vertex] > max_saturated){
                //found vertex with higher saturation, update
                max_saturated = num_selected_neighbors[vertex];
                max_saturated_vertex = vertex;

                saturated_tied_max_degrees.clear();
                saturated_tied_max_degrees.push_back(vertex);
            } else if(num_selected_neighbors[vertex] == max_saturated){
                if(neighbors[vertex].size() > neighbors[max_saturated_vertex].size()){
                    //found vertex with equal saturation but higher degree, update
                    max_saturated_vertex = vertex;

                    saturated_tied_max_degrees.clear();
                    saturated_tied_max_degrees.push_back(vertex);
                } else if(neighbors[vertex].size() == neighbors[max_saturated_vertex].size()){
                    //same saturation level, same degree, add to set of possible candidates
                    saturated_tied_max_degrees.push_back(vertex);
                }
            }
        }
        //select max saturated vertex, at random if there exist multiple
        if(saturated_tied_max_degrees.size() > 1 and random_tiebreaks){
            max_saturated_vertex = random_element(saturated_tied_max_degrees);
        }
        ordering.push_back(max_saturated_vertex);

        int lowest_color = 0;
        bool failed;
        do {
            failed = false;
            for (VertexType neighbor : neighbors[max_saturated_vertex]) {
                if (vertex_color[neighbor] == lowest_color) {
                    ++lowest_color;
                    failed = true;
                    break;
                }
            }
        } while (failed);

        if(lowest_color < int(coloring.size())){
            //we can add chosen vertex to some existing color class
            coloring[lowest_color].insert(max_saturated_vertex);
            vertex_color[max_saturated_vertex] = lowest_color;
        } else if(lowest_color == int(coloring.size())){
            //new color class is created for max_saturated_vertex
            //before that we try to recolor two vertices to not have to use a new color
            bool recolored = try_color_swap(max_saturated_vertex, neighbors, coloring, vertex_color);
            if(not recolored){
                //not able to recolor, add new color class
                coloring.emplace_back();
                coloring[lowest_color].insert(max_saturated_vertex);
                vertex_color[max_saturated_vertex] = lowest_color;
            }
        } else{
            std::cout << "Error: lowest color is too high, something must have gone wrong.\n";
            throw std::runtime_error("Error in Dsatur coloring");
        }

        //new vertex has been selected so update the number of colored neighbors for that vertex's neighbors
        for(VertexType vertex : neighbors[max_saturated_vertex]){
            //possibly check that neighbor has been selected, i.e. it's level is equal to numeric_limits<int>::min()
            if(num_selected_neighbors[vertex] != std::numeric_limits<int>::min()){
                num_selected_neighbors[vertex]++;
            }
        }
        //remove selected vertex from any further consideration
        num_selected_neighbors[max_saturated_vertex] = std::numeric_limits<int>::min();
    }
    return coloring;
}


Permutation Graph::max_connected_degree_ordering(const std::vector<VertexType> &clique) const {
    NeighborList neighbors = get_neighbor_list();
    return max_connected_degree_ordering(neighbors, clique);
}


Permutation Graph::max_connected_degree_ordering(const NeighborList &neighbors, const std::vector<VertexType> &clique) const {
    Permutation ordering;
    ordering.reserve(_ncount);
    //record for an unselected vertex to how many selected vertices it is connected
    std::vector<int> num_selected_neighbors(_ncount);
    std::vector<bool> vertex_is_colored(_ncount, false);

    int nclique_colors = 0;
    if(clique.size() > 2){//ignore trivial or empty clique
        //fix color of the vertices of supplied clique and set saturation levels
        for(VertexType v : clique){
            ordering.push_back(v);
            vertex_is_colored[v] = true;
            //update saturation levels
            for(VertexType neighbor : neighbors[v]){
                num_selected_neighbors[neighbor]++;
            }
            num_selected_neighbors[v] = std::numeric_limits<int>::min();
        }
        nclique_colors = int(clique.size());
    }

    for(int iteration = 0; iteration < _ncount - nclique_colors; iteration++){
        VertexType max_connections_vertex = 0;
        int max_connections = -1;
        for(VertexType vertex = 0; vertex < _ncount; vertex++){
            if(vertex_is_colored[vertex])
                continue;
            if(num_selected_neighbors[vertex] > max_connections){
                max_connections = num_selected_neighbors[vertex];
                max_connections_vertex = vertex;
            } else if(num_selected_neighbors[vertex] == max_connections and
                      (neighbors[vertex].size() > neighbors[max_connections_vertex].size())){
                max_connections_vertex = vertex; //saturation is the same but vertex has a higher degree
            }
        }
        //select max conneceted vertex
        ordering.push_back(max_connections_vertex);
        vertex_is_colored[max_connections_vertex] = true;

        for(VertexType neighbor : neighbors[max_connections_vertex]){
            if(num_selected_neighbors[neighbor] != std::numeric_limits<int>::min()){
                num_selected_neighbors[neighbor]++;
            }
        }
        //remove selected vertex from any further consideration
        num_selected_neighbors[max_connections_vertex] = std::numeric_limits<int>::min();
    }
    return ordering;
}


int Graph::constraint_graph_width() {
    Graph g = *this;
    g.peel_graph(1); //remove isolated vertices
    NeighborList neighbors = g.get_neighbor_list();
    //begin peeling starting from the minimum degree
    auto min_el = std::min_element(neighbors.begin(), neighbors.end(),
                             [](const std::vector<VertexType> &d_1, const std::vector<VertexType> &d_2) {
                                 return d_1.size() < d_2.size();
                             })->size();
    int k = static_cast<int>(min_el);
    while(g.ncount()){
        k++;
        VertexType current_node_count = g.ncount();
        while(true){
            g.peel_graph(k);
            if(current_node_count == g.ncount() or g.ncount() == 0){
                break;
            }
            current_node_count = g.ncount();
        }
    }
    return k;
}

Permutation Graph::constraint_graph_ordering() {
    int width = constraint_graph_width();
    int n = int(_ncount);
    Permutation ordering(n);
    std::set<VertexType> unselected_vertices;
    NeighborList neighbors = get_neighbor_list();
    std::vector<int> degree(_ncount);
    for(VertexType i = 0; i < _ncount; ++i){
        unselected_vertices.insert(unselected_vertices.end(), i);
        degree[i] = int(neighbors[i].size());
    }

    for(int i = n - 1; i >= 0; i--){
        for(auto v_it = unselected_vertices.rbegin(); v_it != unselected_vertices.rend(); v_it++){
            VertexType v = *v_it;
            if(degree[v] <= width + 1){
                ordering[i] = v;
                unselected_vertices.erase(v);
                for(auto neighbor : neighbors[v]){
                    degree[neighbor]--;
                }
                break;
            }
        }
    }
    return ordering;
}


bool Graph::try_color_swap(VertexType max_saturated_vertex, const NeighborList &neighbors, Coloring &coloring,
                           std::vector<int> &vertex_color) {
    for(VertexType j = 0; j < static_cast<int>(coloring.size()); j++){
        for(VertexType k = j + 1; k < static_cast<int>(coloring.size()); k++){
            std::set<VertexType> neighbors_v_intersect_color_j;
            std::set_intersection(neighbors[max_saturated_vertex].begin(),
                                  neighbors[max_saturated_vertex].end(),
                                  coloring[j].begin(), coloring[j].end(),
                                  std::inserter(neighbors_v_intersect_color_j, neighbors_v_intersect_color_j.end()));

            if(neighbors_v_intersect_color_j.size() == 1){
                VertexType u = *neighbors_v_intersect_color_j.begin();
                std::set<VertexType> neighbors_u_intersect_color_k;
                std::set_intersection(neighbors[u].begin(), neighbors[u].end(),
                                      coloring[k].begin(), coloring[k].end(),
                                      std::inserter(neighbors_u_intersect_color_k,
                                                    neighbors_u_intersect_color_k.end()));

                if(neighbors_u_intersect_color_k.empty()){
                    //we can swap colors around, i.e. color max_saturated vertex and u with j and k respectively
                    coloring[j].erase(u);
                    coloring[j].insert(max_saturated_vertex);
                    coloring[k].insert(u);
                    vertex_color[max_saturated_vertex] = int(j);
                    vertex_color[u] = int(k);
                    return true;//was able to swap colors
                }
            }
        }
    }
    return false; //unable to swap two colors
}


Permutation Graph::max_degree_ordering() const {
    NeighborList neighbors = get_neighbor_list();
    return max_degree_ordering(neighbors);
}


Permutation Graph::max_degree_ordering(const NeighborList &neighbors) const {
    Permutation ordering;
    ordering.reserve(_ncount);
    for(VertexType i = 0; i < _ncount; i++){
        ordering.push_back(i);
    }
    std::sort(ordering.begin(), ordering.end(),
              [&neighbors](VertexType v, VertexType w) { return neighbors[v].size() < neighbors[w].size(); });
    std::reverse(ordering.begin(), ordering.end());

    return ordering;
}

void Graph::remove_vertex(VertexType remove) {
    if(remap.empty()){
        //identity mapping
        remap.resize(_ncount);
        for (int v = 0; v < _ncount; ++v) {
            remap[v] = v;
        }
    }
    //update mapping of vertices
    for (int v = 0; v < static_cast<int>(remap.size()); ++v) {
        if(remove == remap[v]){
            remap[v] = UndefVertex; //remove vertex
        }
        else if(remove < remap[v] and remap[v] != UndefVertex){
            remap[v]--; //removed vertex below this one, decrease label of this one
        }
    }
    removed_edges.emplace_back(); //create new vector of edges for this vertex
    //remove and reassign edges
    for(int edge_index = _ecount - 1; edge_index >= 0; edge_index--){
        if(remove == _elist[2 * edge_index] or remove == _elist[2 * edge_index + 1]){
            //remove edges in reverse order
            removed_edges.back().push_back(_elist[2 * edge_index + 1]);
            removed_edges.back().push_back(_elist[2 * edge_index]);
            _elist.erase(std::next(_elist.begin(), 2 * edge_index + 1));
            _elist.erase(std::next(_elist.begin(), 2 * edge_index));
            _ecount--;
        }
        else{
            if(_elist[2 * edge_index + 1] > remove) {
                _elist[2 * edge_index + 1]--;
            }
            if(_elist[2 * edge_index] > remove) {
                _elist[2 * edge_index]--;
            }
        }
    }
    _ncount--;
    removed_vertices.push_back(remove);
    assert(removed_vertices.size() == removed_edges.size());
}


void Graph::remove_vertices(const std::vector<VertexType> &to_remove) {
    assert(std::is_sorted(to_remove.begin(), to_remove.end()));
    //order of removing vertices is important here, by removing later vertices first, we keep the labeling intact
    for(auto v_it = to_remove.rbegin(); v_it != to_remove.rend(); v_it++){
//        std::cout << "removed vertex " << *v_it << "\n";
        remove_vertex(*v_it);
    }
}

//just a tiny helper function to remove vertices all at once
int Graph::num_larger(VertexType v, const std::vector<VertexType>& vec){
    if(vec.empty())
        return 0;
    auto num = std::count_if(vec.begin(), vec.end(), [&v](VertexType set_vertex){return v > set_vertex;});
    return int(num);
}

void Graph::remove_vertices_together(const std::vector<VertexType> &to_remove) {
    assert(std::is_sorted(to_remove.begin(), to_remove.end()));
    if(remap.empty()){
        //identity mapping
        remap.resize(_ncount);
        for (int v = 0; v < _ncount; ++v) {
            remap[v] = v;
        }
    }
    //update mapping of vertices
    for (int v = 0; v < static_cast<int>(remap.size()); ++v) {
        if(std::find(to_remove.begin(), to_remove.end(), remap[v]) != to_remove.end()){
            remap[v] = UndefVertex; //remove vertex
        }
        else if(remap[v] != UndefVertex){
            remap[v] -= num_larger(remap[v], to_remove); //decrease label by number of smaller removed vertices
        }
    }
    //remove and reassign edges
    for (auto rm = to_remove.rbegin(); rm != to_remove.rend(); ++rm) {
        VertexType remove = *rm;
        removed_edges.emplace_back(); //create new vector of edges for this vertex
        for(int edge_index = _ecount - 1; edge_index >= 0; edge_index--) {
            if (remove == _elist[2 * edge_index] or remove == _elist[2 * edge_index + 1]) {
                //remove edges in reverse order
                removed_edges.back().push_back(_elist[2 * edge_index + 1]);
                removed_edges.back().push_back(_elist[2 * edge_index]);
                _elist.erase(std::next(_elist.begin(), 2 * edge_index + 1));
                _elist.erase(std::next(_elist.begin(), 2 * edge_index));
                _ecount--;
            }
            else{
                if(_elist[2 * edge_index + 1] > remove) {
                    _elist[2 * edge_index + 1]--;
                }
                if(_elist[2 * edge_index] > remove) {
                    _elist[2 * edge_index]--;
                }
            }
        }
    }
    _ncount -= static_cast<int>(to_remove.size());
    removed_vertices.insert(removed_vertices.end(), to_remove.rbegin(), to_remove.rend());
    assert(removed_vertices.size() == removed_edges.size());
}

int Graph::peel_graph(int peeling_degree) {
    NeighborList neighbors = get_neighbor_list();
    return peel_graph(peeling_degree, neighbors);
}

int Graph::peel_graph(int peeling_degree, const NeighborList &neighbors) {
    std::vector<VertexType> small_degree_vertices;
    for(VertexType i = 0; i < _ncount; i++){
        if(static_cast<int>(neighbors[i].size()) < peeling_degree){
            small_degree_vertices.push_back(i);
        }
    }
    //vector is already sorted by construction
    remove_vertices_together(small_degree_vertices);
//    remove_vertices(small_degree_vertices);
    return static_cast<int>(small_degree_vertices.size());
}


int Graph::remove_dominated_vertices() {
    NeighborList neighbors = get_neighbor_list();
    return remove_dominated_vertices(neighbors);
}


int Graph::remove_dominated_vertices(const NeighborList &neighbors) {
    std::vector<VertexType> dominated_vertices;
    for(VertexType v = 0; v < _ncount; v++){
        for(VertexType w = v + 1; w < _ncount; w++){
            //vertices are not neighbors
            if(std::find(neighbors[v].begin(), neighbors[v].end(), w) != neighbors[v].end())
                continue;
            //vertices have not been marked as dominated
            if(std::find(dominated_vertices.begin(), dominated_vertices.end(), v) != dominated_vertices.end()
                or std::find(dominated_vertices.begin(), dominated_vertices.end(), w) != dominated_vertices.end()){
                continue;
            }

            std::vector<VertexType> intersection;
            intersection.reserve(std::min(neighbors[v].size(), neighbors[w].size()));
            std::set_intersection(neighbors[v].begin(), neighbors[v].end(),
                                  neighbors[w].begin(), neighbors[w].end(),
                    std::back_inserter(intersection));
            if(intersection == neighbors[v]){
                //neighbors of v included in neighbors of w
                dominated_vertices.push_back(v);
                break;//v has been marked as vertex to be removed
            }
            else if(intersection == neighbors[w]){
                //neighbors of w included in neighbors of v
                dominated_vertices.push_back(w);
            }
        }
    }
    std::sort(dominated_vertices.begin(), dominated_vertices.end());
    remove_vertices_together(dominated_vertices);
//    remove_vertices(dominated_vertices);
    return static_cast<int>(dominated_vertices.size());
}

std::vector<VertexType> Graph::get_remapping() const {
    return remap;
}

std::vector<VertexType> Graph::recover_reductions() {
    assert(not removed_vertices.empty());
    assert(not removed_edges.empty());
    assert(removed_vertices.size() == removed_edges.size());
    assert(not remap.empty());

    //store vertices that were recovered, update them correctly to be of the right vertex labeling
    std::vector<VertexType> recovered_vertices;

    //to insert vertices again, start with last one removed, update vertex mapping
    // and thus edges, then insert removed edges
    //extra iterator to iterate over vectors of removed edges
    auto re_vec = removed_edges.rbegin();
    for (auto rv = removed_vertices.rbegin(); rv != removed_vertices.rend(); ++rv, ++re_vec) {
        //relabel edges to before vertex was removed
        for(int edge_index = _ecount - 1; edge_index >= 0; edge_index--) {
            if (_elist[2 * edge_index + 1] >= *rv) {
                _elist[2 * edge_index + 1]++;
            }
            if (_elist[2 * edge_index] >= *rv) {
                _elist[2 * edge_index]++;
            }
        }
        //insert edges that were removed with removal of that vertex
        for (auto re = (*re_vec).rbegin(); re != (*re_vec).rend(); ++re){
            _elist.push_back(*re);
        }
        _ecount += (static_cast<int>((*re_vec).size()) / 2);
        assert(_ecount == static_cast<int>(_elist.size()/2));

        //update vector of recovered vertices
        for (VertexType &recovered : recovered_vertices) {
            if(recovered >= *rv){
                recovered++;
            }
        }
        recovered_vertices.push_back(*rv);
    }
    _ncount += static_cast<int>(removed_vertices.size());
    assert(_ncount == static_cast<int>(remap.size()));

    removed_vertices.clear();
    removed_edges.clear();
    remap.clear();
    return recovered_vertices;
}


void Graph::write_dimacs(const std::string& out_path) {
    std::ofstream out_file(out_path);
    if(not out_file.is_open()){
        throw std::runtime_error("Unable to open file to write graph to");
    }
    //write header
    out_file << "c Graph in Dimacs format\n";
    out_file << "p edge " << ncount() << " " << ecount() << "\n";
    //write edge list
    for(int i = 0; i < _ecount; i++){
        out_file << "e " <<  _elist[2 * i + 1] + 1 << " " << _elist[2 * i] + 1 << "\n";
    }

    out_file.close();
}


} //namespace Graph



