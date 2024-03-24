#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include "utils.hpp"
#include "Refactor/Refactor.hpp"
#include "reorder.hpp"

template <class T>
void weight_reorder_1D(const T * data_pos, size_t n_nodal, size_t n_coeff, T * nodal_buffer, T * coeff_buffer){
    T * nodal_pos = nodal_buffer;
    T * coeff_pos = coeff_buffer;
    T const * cur_data_pos = data_pos;
    for(int i=0; i<n_coeff; i++){
        *(nodal_pos++) = *(cur_data_pos++);
        *(coeff_pos++) = *(cur_data_pos++);
    }
    *(nodal_pos++) = *(cur_data_pos++);
    if(n_nodal == n_coeff + 2){
        // *nodal_pos = 2*cur_data_pos[0] - nodal_pos[-1];
        // modified process 
        *nodal_pos = *cur_data_pos;
    }
}

/*
    oxoxo       oooxx       oooxx
    xxxxx   (1) xxxxx   (2) oooxx
    oxoxo   =>  oooxx   =>  oooxx
    xxxxx       xxxxx       xxxxx
    oxoxo       oooxx       xxxxx
*/
template <class T>
void weight_reorder_2D(T * data_pos, T * data_buffer, size_t n1, size_t n2, size_t stride){
    size_t n1_nodal = (n1 >> 1) + 1;
    size_t n1_coeff = n1 - n1_nodal;
    size_t n2_nodal = (n2 >> 1) + 1;
    size_t n2_coeff = n2 - n2_nodal;
    T * cur_data_pos = data_pos;
    T * nodal_pos = data_buffer;
    T * coeff_pos = data_buffer + n2_nodal;
    // do reorder (1)
    for(int i=0; i<n1; i++){
        weight_reorder_1D(cur_data_pos, n2_nodal, n2_coeff, nodal_pos, coeff_pos);
        memcpy(cur_data_pos, data_buffer, n2 * sizeof(T));
        cur_data_pos += stride;
    }
    // do reorder (2)
    // TODO: change to online processing for memory saving
    MGARD::switch_rows_2D_by_buffer(data_pos, data_buffer, n1, n2, stride);
}

template <class T>
void weight_reorder_3D(T * data_pos, T * data_buffer, size_t n1, size_t n2, size_t n3, size_t dim0_stride, size_t dim1_stride){
    size_t n1_nodal = (n1 >> 1) + 1;
    size_t n1_coeff = n1 - n1_nodal;
    size_t n2_nodal = (n2 >> 1) + 1;
    size_t n2_coeff = n2 - n2_nodal;
    size_t n3_nodal = (n3 >> 1) + 1;
    size_t n3_coeff = n3 - n3_nodal;
    T * cur_data_pos = data_pos;
    // do 2D reorder
    for(int i=0; i<n1; i++){
        weight_reorder_2D(cur_data_pos, data_buffer, n2, n3, dim1_stride);
        cur_data_pos += dim0_stride;
    }
    cur_data_pos = data_pos;
    // reorder vertically
    for(int j=0; j<n2; j++){
        MGARD::switch_rows_2D_by_buffer(cur_data_pos, data_buffer, n1, n3, dim0_stride);
        cur_data_pos += dim1_stride;
    }
}

using namespace std;
bool negabinary = true;

template <class T, class Refactor>
void evaluate(const vector<T>& data, const vector<uint32_t>& dims, int target_level, int num_bitplanes, Refactor refactor){
    struct timespec start, end;
    int err = 0;
    cout << "Start refactoring" << endl;
    err = clock_gettime(CLOCK_REALTIME, &start);
    refactor.refactor(data.data(), dims, target_level, num_bitplanes);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Refactor time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
}

template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorCollector, class Writer>
void test(string filename, const vector<uint32_t>& dims, int target_level, int num_bitplanes, Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorCollector collector, Writer writer){
    auto refactor = MDR::ComposedRefactor<T, Decomposer, Interleaver, Encoder, Compressor, ErrorCollector, Writer>(decomposer, interleaver, encoder, compressor, collector, writer);
    refactor.negabinary = negabinary;
    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    evaluate(data, dims, target_level, num_bitplanes, refactor);
}

int main(int argc, char ** argv){

    int argv_id = 1;
    string filename = string(argv[argv_id ++]);
    int target_level = atoi(argv[argv_id ++]);
    int num_bitplanes = atoi(argv[argv_id ++]);
    if(num_bitplanes % 2 == 1) {
        num_bitplanes += 1;
        std::cout << "Change to " << num_bitplanes + 1 << " bitplanes for simplicity of negabinary encoding" << std::endl;
    }
    int num_dims = atoi(argv[argv_id ++]);
    vector<uint32_t> dims(num_dims, 0);
    for(int i=0; i<num_dims; i++){
        dims[i] = atoi(argv[argv_id ++]);
    }

    string metadata_file = "refactored_data/metadata.bin";
    vector<string> files;
    for(int i=0; i<=target_level; i++){
        string filename = "refactored_data/level_" + to_string(i) + ".bin";
        files.push_back(filename);
    }
    // using T = float;
    // using T_stream = uint32_t;
    // if(num_bitplanes > 32){
    //     num_bitplanes = 32;
    //     std::cout << "Only less than 32 bitplanes are supported for single-precision floating point" << std::endl;
    // }
    using T = double;
    using T_stream = uint64_t;
    if(num_bitplanes > 64){
        num_bitplanes = 64;
        std::cout << "Only less than 64 bitplanes are supported for double-precision floating point" << std::endl;
    }

    auto decomposer = MDR::MGARDHierarchicalDecomposer<T>();
    auto interleaver = MDR::DirectInterleaver<T>();
    // auto interleaver = MDR::SFCInterleaver<T>();
    // auto interleaver = MDR::BlockedInterleaver<T>();
    // auto encoder = MDR::GroupedBPEncoder<T, T_stream>();
    // auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
    // negabinary = true;
    // auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
    auto encoder = MDR::WeightedPerBitBPEncoder<T, T_stream>();
    negabinary = false;
    auto weight_interleaver = MDR::DirectInterleaver<int>();
    size_t num_elements = 1;
    vector<uint32_t> dims_prev(num_dims);
    for(int i=0; i<dims.size(); i++){
        dims_prev[i] = (dims[i] >> 1) + 1;
        num_elements *= dims[i];
    }
    // int * weights = (int *) malloc(num_elements * sizeof(int));
    // for(int i=0; i<num_elements/2; i++){
    //     weights[i] = 1;
    // }
    auto weights = MGARD::readfile<int>("/Users/xliang/Github/MDR_abs/build/weights_vx.dat", num_elements);
    int max_weight = 0;
    for(int i=0; i<num_elements; i++){
        if(weights[i] > max_weight) max_weight = weights[i];
    }
    std::cout << "max_weight = " << max_weight << std::endl;
    // max_weight = 0;
    int * reordered_weights = (int *) malloc(num_elements * sizeof(int));
    // use reordered_weights as buffer
    weight_reorder_3D(weights.data(), reordered_weights, dims[0], dims[1], dims[2], dims[1]*dims[2], dims[2]);
    weight_interleaver.interleave(weights.data(), dims, dims, dims_prev, reordered_weights);
    // free(weights);
    encoder.set_weights(reordered_weights, max_weight);
    // auto compressor = MDR::DefaultLevelCompressor();
    auto compressor = MDR::AdaptiveLevelCompressor(64);
    // auto compressor = MDR::NullLevelCompressor();
    auto collector = MDR::SquaredErrorCollector<T>();
    auto writer = MDR::ConcatLevelFileWriter(metadata_file, files);
    // auto writer = MDR::HPSSFileWriter(metadata_file, files, 2048, 512 * 1024 * 1024);

    test<T>(filename, dims, target_level, num_bitplanes, decomposer, interleaver, encoder, compressor, collector, writer);
    free(reordered_weights);
    return 0;
}