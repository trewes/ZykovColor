#include "mnts.h"

void MNTS::initializeFromAdjMatrix(const std::vector<Bitset>& adjMatrixInput) {
    Max_Vtx = adjMatrixInput.size(); // Number of vertices
    Iter = 0, Wbest = 0, Wf = 0, len0 = 0, len1 = 0, len_best = 0, len = 0, Titer = 0;

    // Resize and initialize the graph-related data structures
    vectex = Bitset(Max_Vtx);
    funch.resize(Max_Vtx);
    address.resize(Max_Vtx);
    tabuin.resize(Max_Vtx);
    C0.resize(Max_Vtx);
    C1.resize(Max_Vtx);
    TC1.resize(Max_Vtx);
    BC.resize(Max_Vtx);
    FC1.resize(Max_Vtx);
    cruset.resize(2000);
    Tbest = Bitset(Max_Vtx);

    We = std::vector<int>(Max_Vtx, 1); //only unweighted
    Edge = adjMatrixInput; // copy adjacency from input

    // invert edges to obtain complement graph
    for (int x = 0; x < Max_Vtx; x++) {
        Edge[x].flip();
    }
    // remove self loops
    for (int x = 0; x < Max_Vtx; x++) {
        Edge[x].reset(x);
    }
    // Rebuild the adjacency list based on inverted edges
    adjacMatrix.clear();
    adjacMatrix.resize(Max_Vtx);
    for (int x = 0; x < Max_Vtx; x++) {
        for (int y = 0; y < Max_Vtx; y++) {
            if (Edge[x].test(y)) {
                adjacMatrix[x].push_back(y);
            }
        }
    }
}


void MNTS::clearGamma() {
    vectex.reset();
    std::fill(funch.begin(), funch.end(), 0);
    std::fill(address.begin(), address.end(), 0);
    std::fill(tabuin.begin(), tabuin.end(), 0);

    for (int i = 0; i < Max_Vtx; i++) {
        C0[i] = i;
        address[i] = i;
    }
    len0 = Max_Vtx;
    len1 = 0;
    len = 0;
    Wf = 0;
    Wbest = 0;
}

int MNTS::selectC0() {
    if (len0 > 30) {
        return randomInt(len0);
    }
    TC1.clear();
    for (int i = 0; i < len0; i++) {
        int k = C0[i];
        if (tabuin[k] <= Iter) {
            TC1.push_back(i);
        }
    }
    if (TC1.empty()) return -1;

    int index = randomInt(TC1.size());
    return TC1[index];
}

int MNTS::WselectC0() {
    FC1.clear();
    TC1.clear();
    int w1 = 0, w2 = 0;

    for (int i = 0; i < len0; i++) {
        int k = C0[i];
        if (tabuin[k] <= Iter) {
            if (We[k] > w1) {
                w1 = We[k];
                FC1.clear();
                FC1.push_back(i);
            } else if (We[k] == w1) {
                FC1.push_back(i);
            }
        } else {
            if (We[k] > w2) {
                w2 = We[k];
                TC1.clear();
                TC1.push_back(i);
            } else if (We[k] == w2) {
                TC1.push_back(i);
            }
        }
    }

    if (!TC1.empty() && (w2 > w1) && ((w2 + Wf) > Wbest)) {
        return TC1[randomInt(TC1.size())];
    }
    if (!FC1.empty()) {
        return FC1[randomInt(FC1.size())];
    }
    return -1;
}

int MNTS::expand(int SelN) {
    int m = C0[SelN];
    cruset[len++] = m;
    vectex.set(m);
    Wf += We[m];

    len0--;
    int n1 = C0[len0];
    int k1 = address[m];
    C0[k1] = n1;
    address[n1] = k1;

    for (int i = 0; i < adjacMatrix[m].size(); i++) {
        int n = adjacMatrix[m][i];
        funch[n]++;
        if (funch[n] == 1) {
            k1 = address[n];
            len0--;
            n1 = C0[len0];
            C0[k1] = n1;
            address[n1] = k1;

            C1[len1] = n;
            address[n] = len1;
            len1++;
            BC[n] = m;
        } else if (funch[n] == 2) {
            len1--;
            n1 = C1[len1];
            k1 = address[n];
            C1[k1] = n1;
            address[n1] = k1;
        }
    }

    if (Wf > Wbest) {
        Wbest = Wf;
        len_best = len;
        Tbest = vectex;
    }

    return 1;
}



int MNTS::selectC1() {
    TC1.clear();
    for (int i = 0; i < len1; i++) {
        int k = C1[i];
        if (tabuin[k] <= Iter) {
            TC1.push_back(i);
        }
    }
    if (TC1.empty()) return -1;

    int index = randomInt(TC1.size());
    return TC1[index];
}


int MNTS::WselectC1() {
    FC1.clear();
    TC1.clear();
    int w1 = -1000000, w2 = -1000000;
    int l = 0, k = 0;

    for(int i = 0; i < len1; i++ ){
        int m = C1[i];
        int n = BC[m];
        if(vectex.test(n) && Edge[m].test(n)) {
            l++;
        }
        else{
            for(int j = 0; j < len; j++ ){
                k = cruset[j];
                if(Edge[m].test(k)) {
                    break;
                }
            }
            BC[m] = k;
        }
    }

    for(int i = 0; i < len1; i++ )
    {
        int m = C1[i];
        int n = BC[m];
        int wmn = We[ m ] - We[ n ];
        if( tabuin[ m ] <= Iter ){
            if( wmn > w1 ){
                w1 = wmn;
                FC1.clear();
                FC1.push_back(i);
            }
            else if ( wmn == w1 ){
                FC1.push_back(i);
            }
        }
        else{
            if( wmn > w2 ){
                w2 = wmn;
                TC1.clear();
                TC1.push_back(i);
            }
            else if ( wmn == w2 ){
                TC1.push_back(i);
            }
        }
    }

    if (!TC1.empty() && w2 > w1 && (Wf - w2 + Wbest > Waim)) {
        return TC1[randomInt(TC1.size())];
    }
    if (!FC1.empty()) {
        return FC1[randomInt(FC1.size())];
    }
    return -1;
}





int MNTS::plateau(int SelN) {
    int m = C1[SelN];
    int m1 = -1, n = -1, n1 = -1, k1 = -1, ti = -1;

    // Find the first m1 such that Edge[m1][m] == 1
    for (ti = 0; ti < len; ti++) {
        m1 = cruset[ti];
        if (Edge[m1].test(m))
            break;
    }

    // Update Wf by removing m1 and adding m
    Wf = Wf + We[m] - We[m1];

    // The expand process: put m into the current independent set
    vectex.set(m);
    cruset[len++] = m;

    // Remove m from C1
    k1 = address[m];
    len1--;
    n1 = C1[len1];
    C1[k1] = n1;
    address[n1] = k1;

    // Update adjacencies based on the expansion of m
    for (int i = 0; i < adjacMatrix[m].size(); i++) {
        n = adjacMatrix[m][i];
        funch[n]++;
        if (funch[n] == 1 && not vectex.test(n)) {
            // Move n from C0 to C1
            k1 = address[n];
            len0--;
            n1 = C0[len0];
            C0[k1] = n1;
            address[n1] = k1;

            C1[len1] = n;
            address[n] = len1;
            len1++;
            BC[n] = m;
        }
        if (funch[n] == 2) {
            len1--;
            n1 = C1[len1];
            k1 = address[n];
            C1[k1] = n1;
            address[n1] = k1;
        }
    }

    // The backtrack process: remove m1 from the current independent set
    vectex.reset(m1);
    tabuin[m1] = Iter + TABUL + randomInt(len1 + 2);

    // Update the length and the set
    len--;
    cruset[ti] = cruset[len];
    C1[len1] = m1;
    address[m1] = len1;
    len1++;

    // Adjust adjacencies after removing m1
    for (int i = 0; i < adjacMatrix[m1].size(); i++) {
        n = adjacMatrix[m1][i];
        funch[n]--;
        if (funch[n] == 0 && not vectex.test(n)) {
            k1 = address[n];
            len1--;
            n1 = C1[len1];
            C1[k1] = n1;
            address[n1] = k1;

            C0[len0] = n;
            address[n] = len0;
            len0++;
        } else if (funch[n] == 1) {
            C1[len1] = n;
            address[n] = len1;
            len1++;
        }
    }

    // If the current Wf is better than Wbest, update Wbest and store the best solution
    if (Wf > Wbest) {
        Wbest = Wf;
        len_best = len;
        Tbest = vectex;
    }

    return 1;
}




int MNTS::Mumi_Weight() {
    int w1 = 5000000;
    int k = -1;

    // Iterate over cruset and find the element with the smallest weight
    for (int i = 0; i < len; i++) {
        k = cruset[i];
        if (We[k] < w1) {
            w1 = We[k];
            FC1.clear(); // Clear FC1 before adding the new index
            FC1.push_back(i);
        } else if (We[k] == w1) {
            FC1.push_back(i); // If the weight is equal to the current minimum, add to FC1
        }
    }

    // If no valid element was found, return -1
    if (FC1.empty()) {
        return -1;
    }

    // Randomly select an index from the FC1 vector
    k = randomInt(FC1.size());
    k = FC1[k];
    return k;
}




int MNTS::backtrack() {
    // Get the index of the vertex to backtrack
    int ti = Mumi_Weight();
    if (ti == -1) {
        return -1;
    }
    int m1 = cruset[ti];
    Wf -= We[m1];
    vectex.reset(m1);
    tabuin[m1] = Iter + TABUL;
    // Remove m1 from cruset and adjust the C0/C1 and address arrays
    len--;
    cruset[ti] = cruset[len];
    C0[len0] = m1;
    address[m1] = len0;
    len0++;

    // Process adjacent vertices of m1
    for (int i = 0; i < adjacMatrix[m1].size(); i++) {
        int n = adjacMatrix[m1][i];
        funch[n]--;

        if (funch[n] == 0 && not vectex.test(n)) {
            // Move n to C0
            int k1 = address[n];
            len1--;
            int n1 = C1[len1];
            C1[k1] = n1;
            address[n1] = k1;

            C0[len0] = n;
            address[n] = len0;
            len0++;
        } else if (funch[n] == 1) {
            // Move n to C1
            C1[len1] = n;
            address[n] = len1;
            len1++;
        }
    }

    return 1;
}



int MNTS::tabu(int Max_Iter) {
    int am, am1, ww, ww1, ww2, ti, m1;
    Iter = 0;
    clearGamma();

    // First while loop: keep expanding as long as C0 is not exhausted
    while (true){
        am = selectC0();
        if (am != -1){
            int l = expand(am);
            Iter++;
            if (Wbest >= Waim)
                return Wbest;
        }
        else{
            break;
        }
    }

    // Main loop: Continue until we hit Max_Iter or Waim
    while (Iter < Max_Iter){
        am = WselectC0();
        am1 = WselectC1();

        if (am != -1 && am1 != -1){
            ww = We[C0[am]];
            ww1 = We[C1[am1]] - We[BC[C1[am1]]];

            if (ww > ww1)
            {
                int l = expand(am);
                Iter++;
                if (Wbest >= Waim)
                    return Wbest;
            }
            else{
                int l = plateau(am1);
                if (Wbest >= Waim)
                    return Wbest;
                Iter++;
            }
        }
        else if (am != -1 && am1 == -1){
            int l = expand(am);
            if (Wbest >= Waim)
                return Wbest;
            Iter++;
        }
        else if (am == -1 && am1 != -1){
            ti = Mumi_Weight();
            m1 = cruset[ti];
            ww1 = We[C1[am1]] - We[BC[C1[am1]]];
            ww2 = -We[m1];

            if (ww1 > ww2)
            {
                int l = plateau(am1);
                if (Wbest >= Waim)
                    return Wbest;
                Iter++;
            }
            else
            {
                int k = backtrack();
                if (k == -1)
                    return Wbest;
                Iter++;
            }
        }
        else if (am == -1 && am1 == -1){
            int k = backtrack();
            if (k == -1)
                return Wbest;
            Iter++;
        }
    }

    return Wbest;
}


bool MNTS::verify(const Bitset &test) const {
    for (int i = 0; i < Max_Vtx; i++) {
        if (test.test(i) == 1) {
            for (int j = i + 1; j < Max_Vtx; j++) {
                if (test.test(j) == 1 && Edge[i].test(j)) {
                    std::cout << "Hello, there is something wrong!\n";
                    return false;
                }
            }
        }
    }
    return true;
}


void MNTS::outputResults() const {
    std::cout << "Final Results:\n";
    std::cout << "Best Weight: " << Wbest << "\n";
    std::cout << "Solution Vector:\n";
    for (int i = 0; i < Max_Vtx; i++) {
        if (TTbest.test(i)) {
            std::cout << i << " " ;
        }
    }std::cout << "\n";

    if(not verify(TTbest)) {
        std::cout << "Best solution was NOT an independent set!" << "\n";
    }
}


int MNTS::maxTabu() {
    int l, lbest = 0, M_iter = 0;

    for (int i = 0; i < len_time; i++){
        l = tabu(len_improve);
        M_iter += Iter;

        // If the current length is better, update the best length and time
        if (l > lbest){
            lbest = l;
            Titer = M_iter;
            len_W = len_best;
            TTbest = Tbest;
            assert(verify(TTbest));
        }

        // If the desired length is achieved, return the best found
        if (l >= Waim){
            return lbest;
        }
    }

    return lbest;
}


Bitset MNTS::runSearch() {
    Wbest = maxTabu();
    // outputResults();
    return TTbest;
}
























