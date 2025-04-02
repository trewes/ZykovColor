# ZykovColor


This repository contains the code for the exact graph coloring algorithm **ZykovColor**, 
presented in "A Customized SAT-based Solver for Graph Coloring".

The focus is a custom external propagator using the IPASIR-UP interface [[1]](#1)
to communicate information with the underlying SAT solver CaDiCal [[2]](#2).
It works on top of a Zykov-based SAT encoding for the graph coloring problem
and handles the transitivity constraints, can prune parts of the search with lower bounds,
and uses further techniques detailed in the paper to yield a competitive algorithm.
The presented propagator extends the ideas of [[3]](#3) and uses SAT instead of a CSP solver, 
which also allows solving the decision problems incrementally and re-using information and learned clauses.

## Dependencies and Installation


We used version 2.1.2 of [CaDiCal](https://github.com/arminbiere/cadical), 
but any version with the same IPASIR-UP interface should work.
[Open-WBO](https://github.com/sat-group/open-wbo) is used
to build the incremental cardinality constraints with the Totalizer Encoding.
This is not strictly necessary for ZykovColor, 
but we also support the option to explicitly add at-most-$k$-constraints.
[Boost](https://www.boost.org/) and ZLIB are also required,
and the path to a binary of [CliSAT](https://github.com/psanse/CliSAT) is needed
to compute an initial clique in preprocessing.


To create and use the binary ``ZykovColor`` it should be as simple as
```
mkdir build && cd build
cmake .. -DBOOST_ROOT=/path/to/boost -DOPENWBO_ROOT_DIR=/path/to/open-wbo -DCADICAL_ROOT_DIR=/path/to/cadical -DCLISAT_BINARY_PATH=/path/to/clisat/binary
make
```
Additionally add ``-DCMAKE_BUILD_TYPE=Release`` to compile in optimized mode.

## Usage

```
./IncSatGC inputfile [options] --configuration 
```
Where ``--configuration`` is one of the following:
- ``--zykov-color``
- ``--assignment``
- ``--partial-order``

They run either ZykovColor, the Assignment encoding, or the Partial Order encoding 
in their default configuration as described in the paper.

### Reference

<a id="1">[1]</a>
Fazekas, Katalin, et al. "IPASIR-UP: user propagators for CDCL."
26th International Conference on Theory and Applications of Satisfiability Testing (SAT 2023).
Schloss Dagstuhl-Leibniz-Zentrum für Informatik, 2023.

<a id="1">[2]</a>
Biere, Armin, et al. "CaDiCaL 2.0."
International Conference on Computer Aided Verification.
Cham: Springer Nature Switzerland, 2024.

<a id="1">[3]</a>
Hebrard, Emmanuel, and George Katsirelos. "Constraint and satisfiability reasoning for graph coloring." 
Journal of Artificial Intelligence Research 69 (2020): 33-65.
