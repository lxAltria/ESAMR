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
void evaluate(const vector<T>& data, double tolerance, Reconstructor reconstructor){
    struct timespec start, end;
    int err = 0;

    cout << "Start reconstruction" << endl;
    err = clock_gettime(CLOCK_REALTIME, &start);
    auto reconstructed_data = reconstructor.reconstruct(tolerance);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Reconstruct time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    // TODO: add full resolution check
    MGARD::print_statistics(data.data(), reconstructed_data, data.size());
}

template <class T, class Decomposer, class Interleaver, class Encoder, class SizeInterpreter, class Retriever>
void test(string filename, double tolerance, Decomposer decomposer, Interleaver interleaver, Encoder encoder, SizeInterpreter interpreter, Retriever retriever){
    auto reconstructor = MDR::ComposedReconstructor<T, Decomposer, Interleaver, Encoder, SizeInterpreter, Retriever>(decomposer, interleaver, encoder, interpreter, retriever);
    cout << "loading metadata" << endl;
    reconstructor.load_metadata();

    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    evaluate(data, tolerance, reconstructor);
}

int main(int argc, char ** argv){

    int argv_id = 1;
    string filename = string(argv[argv_id ++]);
    double tolerance = atof(argv[argv_id ++]);

    string metadata_file = "metadata.bin";
    int num_levels = 0;
    int num_dims = 0;
    {
        size_t num_bytes = 0;
        auto metadata = MGARD::readfile<uint8_t>(metadata_file.c_str(), num_bytes);
        assert(num_bytes > num_dims * sizeof(uint32_t) + 2);
        // metadata interpreter, otherwise information needs to be provided
        num_dims = metadata[0];
        num_levels = metadata[num_dims * sizeof(uint32_t) + 1];
        cout << "number of dimension = " << num_dims << ", number of levels = " << num_levels << endl;
    }
    vector<string> files;
    for(int i=0; i<num_levels; i++){
        string filename = "refactored_data/level_" + to_string(i) + ".bin";
        files.push_back(filename);
    }

    using T = float;
    using T_stream = uint32_t;
    auto decomposer = MDR::MGARDOrthoganalDecomposer<T>();
    auto interleaver = MDR::DirectInterleaver<T>();
    auto encoder = MDR::GroupedBPEncoder<T, T_stream>();
    auto estimator = MDR::SNormErrorEstimator<T>(num_dims, num_levels - 1, 0);
    auto interpreter = MDR::InorderSizeInterpreter<MDR::SNormErrorEstimator<T>>(estimator);
    auto retriever = MDR::FileRetriever(metadata_file, files);

    test<T>(filename, tolerance, decomposer, interleaver, encoder, interpreter, retriever);
    return 0;
}