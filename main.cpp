#include <iostream>
#include "IncSatGC.h"

int main(int argc, char *argv[]) {
    Options opt(argc, argv);
    IncSatGC instance(opt.filepath.c_str(), opt);
    int chromatic_number = instance.run();
    if(opt.strategy == Options::SingleK) {
        std::cout << "IncSatGC computed " << (chromatic_number ? "SAT" : "UNSAT") << " for k = " << opt.specific_num_colors.value() << "\n";
    }
    else{
        std::cout << "IncSatGC computed chromatic number of " << chromatic_number << "\n";
    }
    return 0;
}