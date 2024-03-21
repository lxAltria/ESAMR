#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <bitset>
#include <numeric>
#include "utils.hpp"
#include "Synthesizer4GE.hpp"

using namespace std;

int main(int argc, char** argv){

    using T = double;
    int arg_pos = 1;
    int id = atoi(argv[arg_pos++]);

    MDR::test_singleZone<T>(id);

    return 0;

}