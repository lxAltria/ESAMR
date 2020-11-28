#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include "utils.hpp"
#include "Refactor/Refactor.hpp"

using namespace std;

template <class T, class Refactor>
void evaluate(const vector<T>& data, const vector<uint32_t>& dims, int target_level, int num_bitplanes, Refactor refactor){
    struct timespec start, end;
    int err = 0;

    uint32_t refactored_size = 0;
    cout << "Start refactoring" << endl;
    err = clock_gettime(CLOCK_REALTIME, &start);
    auto refactored_data = refactor.refactor(data.data(), dims, target_level, num_bitplanes, refactored_size);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Refactor time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    cout << "Refactored data size = " << refactored_size << endl;
    MGARD::writefile("refactored_data.bin", refactored_data, refactored_size);
    free(refactored_data);

    uint32_t metadata_size = 0;
    auto metadata = refactor.dump_metadata(metadata_size);
    MGARD::writefile("metadata.bin", metadata, metadata_size);
    free(metadata);

}

template <class T, class Decomposer, class Interleaver, class Encoder, class ErrorCollector, class Reorganizer>
void test(string filename, const vector<uint32_t>& dims, int target_level, int num_bitplanes, Decomposer decomposer, Interleaver interleaver, Encoder encoder, ErrorCollector collector, Reorganizer reorganizer){
    auto refactor = MDR::ComposedRefactor<T, Decomposer, Interleaver, Encoder, ErrorCollector, Reorganizer>(decomposer, interleaver, encoder, collector, reorganizer);
    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    evaluate(data, dims, target_level, num_bitplanes, refactor);
}

int main(int argc, char ** argv){

    int argv_id = 1;
    string filename = string(argv[argv_id ++]);
    int target_level = atoi(argv[argv_id ++]);
    int num_bitplanes = atoi(argv[argv_id ++]);
    int num_dims = atoi(argv[argv_id ++]);
    vector<uint32_t> dims(num_dims, 0);
    for(int i=0; i<num_dims; i++){
        dims[i] = atoi(argv[argv_id ++]);
    }

    using T = float;
    using T_stream = uint32_t;
    auto decomposer = MDR::MGARDOrthoganalDecomposer<T>();
    auto interleaver = MDR::DirectInterleaver<T>();
    auto encoder = MDR::GroupedBPEncoder<T, T_stream>();
    auto collector = MDR::SNormErrorCollector<T>();
    // auto collector = MDR::MaxErrorCollector<T>();
    auto reorganizer = MDR::InOrderReorganizer();

    test<T>(filename, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, collector, reorganizer);
    return 0;
}