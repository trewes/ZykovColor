#ifndef FRACTIONALBOUND_H
#define FRACTIONALBOUND_H

#include <unordered_map>
#include <cfenv>
#include <iomanip>

#include "GraphMatrix.h"

#include "boost/dynamic_bitset.hpp"

extern "C" {
#include "color.h"
#include "color_parms.h"
#include "color_private.h"
}

using Bitset = boost::dynamic_bitset<>;

double fractional_chromatic_number_exactcolors(const std::vector<Bitset>& adj);

#endif //FRACTIONALBOUND_H