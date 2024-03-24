#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include "utils.hpp"
#include "Reconstructor/Reconstructor.hpp"
#include "reorder.hpp"
// #include "evaluate.hpp"

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

template <class T, class Reconstructor>
void evaluate(const vector<T>& data, const vector<double>& tolerance, Reconstructor reconstructor){
    struct timespec start, end;
    int err = 0;
    // auto a1 = compute_average(data.data(), dims[0], dims[1], dims[2], 3);
    // auto a12 = compute_average(data.data(), dims[0], dims[1], dims[2], 5);
    for(int i=0; i<tolerance.size(); i++){
        cout << "Start reconstruction" << endl;
        err = clock_gettime(CLOCK_REALTIME, &start);
        auto reconstructed_data = reconstructor.progressive_reconstruct(tolerance[i], -1);
        err = clock_gettime(CLOCK_REALTIME, &end);
        cout << "Reconstruct time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
        auto dims = reconstructor.get_dimensions();
        MGARD::print_statistics(data.data(), reconstructed_data, data.size());
        // COMP_UTILS::evaluate_gradients(data.data(), reconstructed_data, dims[0], dims[1], dims[2]);
        // COMP_UTILS::evaluate_average(data.data(), reconstructed_data, dims[0], dims[1], dims[2], 0);
    }
}

template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorEstimator, class SizeInterpreter, class Retriever>
void test(string filename, const vector<double>& tolerance, Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorEstimator estimator, SizeInterpreter interpreter, Retriever retriever){
    auto reconstructor = MDR::ComposedReconstructor<T, Decomposer, Interleaver, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever>(decomposer, interleaver, encoder, compressor, interpreter, retriever);
    cout << "loading metadata" << endl;
    reconstructor.load_metadata();

    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    std::cout << "read file done: #element = " << num_elements << std::endl;
    fflush(stdout);
    evaluate(data, tolerance, reconstructor);
}

int main(int argc, char ** argv){

    int argv_id = 1;
    string filename = string(argv[argv_id ++]);
    int error_mode = atoi(argv[argv_id++]);
    int num_tolerance = atoi(argv[argv_id ++]);
    vector<double> tolerance(num_tolerance, 0);
    for(int i=0; i<num_tolerance; i++){
        tolerance[i] = atof(argv[argv_id ++]);  
    }
    double s = atof(argv[argv_id ++]);

    string metadata_file = "refactored_data/metadata.bin";
    int num_levels = 0;
    int num_dims = 0;
    vector<uint32_t> dims;
    {
        // metadata interpreter, otherwise information needs to be provided
        size_t num_bytes = 0;
        auto metadata = MGARD::readfile<uint8_t>(metadata_file.c_str(), num_bytes);
        assert(num_bytes > num_dims * sizeof(uint32_t) + 2);
        num_dims = metadata[0];
        num_levels = metadata[num_dims * sizeof(uint32_t) + 1];
        cout << "number of dimension = " << num_dims << ", number of levels = " << num_levels << endl;
        uint8_t * metadata_pos = &metadata[1];
        for(int i=0; i<num_dims; i++){
            dims.push_back(*reinterpret_cast<uint32_t*>(metadata_pos));
            metadata_pos += sizeof(uint32_t);
            cout << dims[i] << " ";
        }
        cout << endl;
        fflush(stdout);
    }
    vector<string> files;
    for(int i=0; i<num_levels; i++){
        string filename = "refactored_data/level_" + to_string(i) + ".bin";
        files.push_back(filename);
    }

    using T = double;
    using T_stream = uint64_t;
    auto decomposer = MDR::MGARDHierarchicalDecomposer<T>();
    auto interleaver = MDR::DirectInterleaver<T>();
    // auto interleaver = MDR::SFCInterleaver<T>();
    // auto interleaver = MDR::BlockedInterleaver<T>();
    // auto encoder = MDR::GroupedBPEncoder<T, T_stream>();
    // auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
    // auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
    auto encoder = MDR::WeightedPerBitBPEncoder<T, T_stream>();
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
    std::cout << dims.size() << " " << dims_prev.size() << std::endl;
    int * reordered_weights = (int *) malloc(num_elements * sizeof(int));
    weight_reorder_3D(weights.data(), reordered_weights, dims[0], dims[1], dims[2], dims[1]*dims[2], dims[2]);
    weight_interleaver.interleave(weights.data(), dims, dims, dims_prev, reordered_weights);
    // free(weights);
    encoder.set_weights(reordered_weights, max_weight);
    // auto compressor = MDR::DefaultLevelCompressor();
    auto compressor = MDR::AdaptiveLevelCompressor(64);
    // auto compressor = MDR::NullLevelCompressor();
    auto retriever = MDR::ConcatLevelFileRetriever(metadata_file, files);
    switch(error_mode){
        case 1:{
            auto estimator = MDR::SNormErrorEstimator<T>(num_dims, num_levels - 1, s);
            // auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::SNormErrorEstimator<T>>(estimator);
            auto interpreter = MDR::NegaBinaryGreedyBasedSizeInterpreter<MDR::SNormErrorEstimator<T>>(estimator);
            // auto interpreter = MDR::RoundRobinSizeInterpreter<MDR::SNormErrorEstimator<T>>(estimator);
            // auto interpreter = MDR::InorderSizeInterpreter<MDR::SNormErrorEstimator<T>>(estimator);
            // auto estimator = MDR::L2ErrorEstimator_HB<T>(num_dims, num_levels - 1);
            // auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::L2ErrorEstimator_HB<T>>(estimator);
            test<T>(filename, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);            
            break;
        }
        default:{
            // auto estimator = MDR::MaxErrorEstimatorOB<T>(num_dims);
            // auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorOB<T>>(estimator);
            // auto interpreter = MDR::RoundRobinSizeInterpreter<MDR::MaxErrorEstimatorOB<T>>(estimator);
            // auto interpreter = MDR::InorderSizeInterpreter<MDR::MaxErrorEstimatorOB<T>>(estimator);
            auto estimator = MDR::MaxErrorEstimatorHB<T>();
            auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>(estimator);
            test<T>(filename, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
        }
    }    
    free(reordered_weights);
    return 0;
}