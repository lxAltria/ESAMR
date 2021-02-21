#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include "utils.hpp"
#include "Reconstructor/Reconstructor.hpp"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

using namespace std;

template <class T>
void posix_read(std::string filename, T * data, size_t num_elements){
  int fd = open(filename.c_str(), O_RDONLY);
  read(fd, data, num_elements * sizeof(T));
  close(fd);
}

template <class T, class Reconstructor>
void evaluate(const vector<T>& data, int num_tolerance, Reconstructor reconstructor){
    struct timespec start, end;
    int err = 0;

    for(int i=0; i<num_tolerance; i++){
        string fragment_file = "refactored_data/fragment_" + to_string(i) + ".bin";
        struct stat st;
        stat(fragment_file.c_str(), &st);
        uint32_t fragment_size = st.st_size;
        uint8_t * fragment_data = (uint8_t *) malloc(fragment_size);
        // cout << "Start reconstruction" << endl;
        err = clock_gettime(CLOCK_REALTIME, &start);
        // read in next fragment
        posix_read(fragment_file, fragment_data, fragment_size);
        // cout << "Read fragment done" << endl;
        auto reconstructed_data = reconstructor.progressive_reconstruct(fragment_data, fragment_size);
        err = clock_gettime(CLOCK_REALTIME, &end);
        cout << "Reconstruct time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
        // TODO: add full resolution check
        free(fragment_data);
        MGARD::print_statistics(data.data(), reconstructed_data, data.size());        
    }
}

template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorEstimator, class SizeInterpreter, class Retriever>
void test(string filename, int num_tolerance, Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorEstimator estimator, SizeInterpreter interpreter, Retriever retriever){
    auto reconstructor = MDR::ComposedReconstructor<T, Decomposer, Interleaver, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever>(decomposer, interleaver, encoder, compressor, interpreter, retriever);
    // cout << "loading metadata" << endl;
    string metadata_file = "refactored_data/metadata.bin";
    struct stat st;
    stat(metadata_file.c_str(), &st);
    uint32_t metadata_size = st.st_size;
    uint8_t * metadata = (uint8_t *) malloc(metadata_size);
    posix_read(metadata_file, metadata, metadata_size);
    reconstructor.load_metadata(metadata);
    free(metadata);

    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    evaluate(data, num_tolerance, reconstructor);
}

int main(int argc, char ** argv){

    string filename = string(argv[1]);
    int num_tolerance = atoi(argv[2]);
    string metadata_file = "refactored_data/metadata.bin";
    int num_levels = 0;
    int num_dims = 0;
    {
        // metadata interpreter, otherwise information needs to be provided
        size_t num_bytes = 0;
        auto metadata = MGARD::readfile<uint8_t>(metadata_file.c_str(), num_bytes);
        assert(num_bytes > num_dims * sizeof(uint32_t) + 2);
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
    // auto decomposer = MDR::MGARDOrthoganalDecomposer<T>();
    auto decomposer = MDR::MGARDHierarchicalDecomposer<T>();
    auto interleaver = MDR::DirectInterleaver<T>();
    // auto interleaver = MDR::SFCInterleaver<T>();
    // auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
    auto encoder = MDR::GroupedBPEncoder<T, T_stream>();
    auto compressor = MDR::DefaultLevelCompressor();
    // auto compressor = MDR::NullLevelCompressor();
    auto retriever = MDR::ConcatLevelFileRetriever(metadata_file, files);
    // auto estimator = MDR::MaxErrorEstimatorOB<T>(num_dims);
    // auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorOB<T>>(estimator);
    auto estimator = MDR::MaxErrorEstimatorHB<T>();
    auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>(estimator);
    test<T>(filename, num_tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
    return 0;
}