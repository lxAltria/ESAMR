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

template <class T>
vector<T> compute_average(const T * data, size_t n1, size_t n2, size_t n3, int block_size){
    size_t dim0_offset = n2 * n3;
    size_t dim1_offset = n3;
    size_t num_block_1 = (n1 - 1) / block_size + 1;
    size_t num_block_2 = (n2 - 1) / block_size + 1;
    size_t num_block_3 = (n3 - 1) / block_size + 1;
    vector<T> aggregated = vector<T>();
    size_t index = 0;
    const T * data_x_pos = data;
    for(int i=0; i<num_block_1; i++){
        int size_1 = (i == num_block_1 - 1) ? n1 - i * block_size : block_size;
        const T * data_y_pos = data_x_pos;
        for(int j=0; j<num_block_2; j++){
            int size_2 = (j == num_block_2 - 1) ? n2 - j * block_size : block_size;
            const T * data_z_pos = data_y_pos;
            for(int k=0; k<num_block_3; k++){
                int size_3 = (k == num_block_3 - 1) ? n3 - k * block_size : block_size;
                const T * cur_data_pos = data_z_pos;
                int n_block_elements = size_1 * size_2 * size_3;
                double sum = 0;
                for(int ii=0; ii<size_1; ii++){
                    for(int jj=0; jj<size_2; jj++){
                        for(int kk=0; kk<size_3; kk++){
                            sum += *cur_data_pos;
                            cur_data_pos ++;
                        }
                        cur_data_pos += dim1_offset - size_3;
                    }
                    cur_data_pos += dim0_offset - size_2 * dim1_offset;
                }
                aggregated.push_back(sum / n_block_elements);
                data_z_pos += size_3;
            }
            data_y_pos += dim1_offset * size_2;
        }
        data_x_pos += dim0_offset * size_1;
    }    
    return aggregated;
}

template <class T>
void compare(const vector<T>& v1, const vector<T>& v2){
    assert(v1.size() == v2.size());
    T max_err = 0;
    T max_v = 0;
    double L2_err = 0;
    for(int i=0; i<v1.size(); i++){
        if(fabs(v1[i] - v2[i]) > max_err) max_err = fabs(v1[i] - v2[i]);
        if(fabs(v1[i]) > max_v) max_v = fabs(v1[i]);
        L2_err += (v1[i] - v2[i]) * (v1[i] - v2[i]);
    }
    cout << "Relative L^infty error = " << max_err / max_v << ", ";
    cout << "L^2 error = " << sqrt(L2_err/v1.size()) << endl;
}

template <class T, class Reconstructor>
void evaluate(const vector<T>& data, const vector<double>& tolerance, Reconstructor reconstructor){
    struct timespec start, end;
    int err = 0;
    auto dims = reconstructor.get_dimensions();
    // auto a1 = compute_average(data.data(), dims[0], dims[1], dims[2], 3);
    // auto a12 = compute_average(data.data(), dims[0], dims[1], dims[2], 5);
    for(int i=0; i<tolerance.size(); i++){
        cout << "Start reconstruction" << endl;
        err = clock_gettime(CLOCK_REALTIME, &start);
        auto reconstructed_data = reconstructor.progressive_reconstruct(tolerance[i]);
        err = clock_gettime(CLOCK_REALTIME, &end);
        cout << "Reconstruct time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
        // TODO: add full resolution check
        MGARD::print_statistics(data.data(), reconstructed_data, data.size());
        // auto a2 = compute_average(reconstructed_data, dims[0], dims[1], dims[2], 3);
        // cout << "Average: ";
        // compare(a1, a2);
        // auto a22 = compute_average(reconstructed_data, dims[0], dims[1], dims[2], 5);
        // cout << "Average2: ";
        // compare(a12, a22);
    }
}

template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorEstimator, class SizeInterpreter, class Retriever>
void test(string filename, const vector<double>& tolerance, Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorEstimator estimator, SizeInterpreter interpreter, Retriever retriever){
    auto reconstructor = MDR::ComposedReconstructor<T, Decomposer, Interleaver, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever>(decomposer, interleaver, encoder, compressor, interpreter, retriever);
    cout << "loading metadata" << endl;
    reconstructor.load_metadata();

    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
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
    auto decomposer = MDR::MGARDOrthoganalDecomposer<T>();
    // auto decomposer = MDR::MGARDHierarchicalDecomposer<T>();
    // auto interleaver = MDR::DirectInterleaver<T>();
    // auto interleaver = MDR::SFCInterleaver<T>();
    auto interleaver = MDR::BlockedInterleaver<T>();
    // auto encoder = MDR::GroupedBPEncoder<T, T_stream>();
    auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
    // auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
    // auto compressor = MDR::DefaultLevelCompressor();
    auto compressor = MDR::AdaptiveLevelCompressor();
    // auto compressor = MDR::NullLevelCompressor();
    auto retriever = MDR::ConcatLevelFileRetriever(metadata_file, files);
    switch(error_mode){
        case 1:{
            auto estimator = MDR::SNormErrorEstimator<T>(num_dims, num_levels - 1, s);
            // auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::SNormErrorEstimator<T>>(estimator);
            auto interpreter = MDR::NegaBinaryGreedyBasedSizeInterpreter<MDR::SNormErrorEstimator<T>>(estimator);
            // auto estimator = MDR::L2ErrorEstimator_HB<T>(num_dims, num_levels - 1);
            // auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::L2ErrorEstimator_HB<T>>(estimator);
            test<T>(filename, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);            
            break;
        }
        default:{
            auto estimator = MDR::MaxErrorEstimatorOB<T>(num_dims);
            auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorOB<T>>(estimator);
            // auto estimator = MDR::MaxErrorEstimatorHB<T>();
            // auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>(estimator);
            test<T>(filename, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
        }
    }    
    return 0;
}