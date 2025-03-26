#ifndef MNTS_CPP_H
#define MNTS_CPP_H

/*  Multi-neighborhood tabu search for the maximum weight clique problem
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*
* This program demonstrates the use of the heuristic algorithm for solving
* the problem of the maximum weight clique problem. For more information about this algorithm,
* please visit the website http://www.info.univ-angers.fr/pub/hao/mnts_clique.html or email to: wu@info-univ.angers.fr.
*/

/*
* The source code was ported from the original c code provided at https://leria-info.univ-angers.fr/~jinkao.hao/clique.html
* to c++ and using the dynamics_bitset data structure where possible
*/

#include <vector>
#include <algorithm>
#include <random>
#include <cstring>
#include <iostream>
#include <fstream>
#include "boost/dynamic_bitset.hpp"
using Bitset = boost::dynamic_bitset<>;

class MNTS {
    int Max_Vtx;
    int Iter;
    int Wbest;
    int Wf;
    int len0, len1, len, len_best;
    int TABUL = 7;
    int Waim;
    int len_time, len_improve;
    int Titer;
    int len_W;

    std::vector<int> We;
    Bitset vectex;
    std::vector<int> funch;
    std::vector<int> address;
    std::vector<int> tabuin;
    std::vector<int> C0, C1, cruset;
    std::vector<int> TC1, FC1;
    std::vector<int> BC;
    std::vector<std::vector<int>> adjacMatrix;
    std::vector<Bitset> Edge;
    Bitset Tbest, TTbest;
    std::vector<Bitset> best_list; //possibly track multiple cliques that achieve the aim, currently not done

    std::mt19937 rng;
    std::uniform_int_distribution<int> dist;

    int randomInt(int range) {
        return dist(rng) % range;
    }

    void initializeFromAdjMatrix(const std::vector<Bitset>& adjMatrixInput);
    void clearGamma();
    int selectC0();
    int WselectC0();
    int expand(int SelN);
    int selectC1();
    int WselectC1();
    int plateau(int SelN);
    int Mumi_Weight();
    int backtrack();
    int tabu(int Max_Iter);

    [[nodiscard]] bool verify(const Bitset &test) const;
    void outputResults() const;
    int maxTabu();

public:

    MNTS(const std::vector<Bitset>& adjMatrixInput, const int Waim,
        const int mnts_length, const int lenImprove, int random_seed = 12345)
        : Waim(Waim), len_improve(lenImprove) {
        initializeFromAdjMatrix(adjMatrixInput);
        // rng.seed(static_cast<unsigned>(std::time(0)));
        rng.seed(random_seed);
        dist = std::uniform_int_distribution<int>();
        len_time = static_cast<int>(mnts_length / len_improve) + 1;
    }
    Bitset runSearch();
};

#endif //MNTS_CPP_H
