#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include "utils.hpp"
#include "Reconstructor/Reconstructor.hpp"

using namespace std;

template <class T, class Reconstructor>
void evaluate(const vector<T>& data, string metadata_file, string data_file, double tolerance, Reconstructor reconstructor){
    struct timespec start, end;
    int err = 0;

    size_t num_bytes = 0;
    auto metadata = MGARD::readfile<uint8_t>(metadata_file.c_str(), num_bytes);
    reconstructor.load_metadata(metadata.data());
    auto refactored_data = MGARD::readfile<uint8_t>(data_file.c_str(), num_bytes);
    cout << "Start reconstruction" << endl;
    err = clock_gettime(CLOCK_REALTIME, &start);
    auto reconstructed_data = reconstructor.reconstruct(refactored_data.data(), tolerance);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Reconstruct time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    // TODO: add full resolution check
    MGARD::print_statistics(data.data(), reconstructed_data, data.size());
}

template <class T, class Decomposer, class Interleaver, class Encoder, class Retriever>
void test(string filename, double tolerance, Decomposer decomposer, Interleaver interleaver, Encoder encoder, Retriever retriever){
    auto reconstructor = MDR::ComposedReconstructor<T, Decomposer, Interleaver, Encoder, Retriever>(decomposer, interleaver, encoder, retriever);
    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    evaluate(data, "metadata.bin", "refactored_data.bin", tolerance, reconstructor);
}

int main(int argc, char ** argv){

    int argv_id = 1;
    string filename = string(argv[argv_id ++]);
    double tolerance = atof(argv[argv_id ++]);

    using T = float;
    using T_stream = uint32_t;
    auto decomposer = MDR::MGARDOrthoganalDecomposer<T>();
    auto interleaver = MDR::DirectInterleaver<T>();
    auto encoder = MDR::GroupedBPEncoder<T, T_stream>();
    // TODO: automate dimensions
    auto estimator = MDR::SNormErrorEstimator<T>(3, 3, 0);
    auto retriever = MDR::GeneralRetriever<MDR::SNormErrorEstimator<T>>(estimator);

    test<T>(filename, tolerance, decomposer, interleaver, encoder, retriever);
    return 0;
}