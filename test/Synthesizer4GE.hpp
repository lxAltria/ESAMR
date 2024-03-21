#ifndef _MDR_GE_SYNTHESIZER_HPP
#define _MDR_GE_SYNTHESIZER_HPP

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cassert>
#include <vector>
#include <cmath>
#include <numeric>
#include "Reconstructor/Reconstructor.hpp"

const std::string data_file_prefix = "/Users/xliang/Dataset/GE/data/";
const std::string rdata_file_prefix = "/Users/xliang/Dataset/GE/refactor/";
const std::vector<std::string> varlist = {"Pressure", "Density", "VelocityX", "VelocityY", "VelocityZ"};
const int n_vars = 5;

namespace MDR {

template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorEstimator, class SizeInterpreter, class Retriever>
MDR::ComposedReconstructor<T, Decomposer, Interleaver, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever> generateReconstructor(Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorEstimator estimator, SizeInterpreter interpreter, Retriever retriever){
    auto reconstructor = MDR::ComposedReconstructor<T, Decomposer, Interleaver, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever>(decomposer, interleaver, encoder, compressor, interpreter, retriever);
    return reconstructor;
}

template<class Type>
void test_singleZone(int zone_id){

    std::string id_str = std::to_string(zone_id);

    size_t num_elements = 0;
    auto pressure_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_Pressure.dat").c_str(), num_elements);
    auto density_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_Density.dat").c_str(), num_elements);
    auto velocityX_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_VelocityX.dat").c_str(), num_elements);
    auto velocityY_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_VelocityY.dat").c_str(), num_elements);
    auto velocityZ_vec = MGARD::readfile<Type>((data_file_prefix + "zone_" + id_str + "_VelocityZ.dat").c_str(), num_elements);

    std::vector<std::vector<Type>> vars_vec = {pressure_vec, density_vec, velocityX_vec, velocityY_vec, velocityZ_vec};

    std::vector<MDR::ComposedReconstructor<Type, MGARDHierarchicalDecomposer<Type>, DirectInterleaver<Type>, PerBitBPEncoder<Type, uint32_t>, AdaptiveLevelCompressor, SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<Type>>, MaxErrorEstimatorHB<Type>, ConcatLevelFileRetriever>> reconstructors;

    for(int i=0; i<n_vars; i++){
        // std::string rdir_prefix = rdata_file_prefix + "zone_" + id_str + "_" + varlist[i];
        std::string rdir_prefix = rdata_file_prefix + varlist[i];
        std::string metadata_file = rdir_prefix + "_refactored_data/metadata.bin";
        std::vector<std::string> files;
        int num_levels = 5;
        for(int i=0; i<num_levels; i++){
            std::string filename = rdir_prefix + "_refactored_data/level_" + std::to_string(i) + ".bin";
            files.push_back(filename);
        }
        auto decomposer = MGARDHierarchicalDecomposer<Type>();
        auto interleaver = DirectInterleaver<Type>();
        auto encoder = PerBitBPEncoder<Type, uint32_t>();
        auto compressor = AdaptiveLevelCompressor(64);
        auto estimator = MaxErrorEstimatorHB<Type>();
        auto interpreter = SignExcludeGreedyBasedSizeInterpreter<MaxErrorEstimatorHB<Type>>(estimator);
        auto retriever = ConcatLevelFileRetriever(metadata_file, files);
        // auto reconstructor = ComposedReconstructor(decomposer, interleaver, encoder, compressor, interpreter, retriever);  // correction
        reconstructors.push_back(generateReconstructor<Type>(decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever));
        reconstructors.back().load_metadata();
        // size_t num_bytes;
        // auto metadata = MGARD::readfile<uint8_t>(metadata_file.c_str(), num_bytes);
        // auto num_dims = metadata[0];
        // num_levels = metadata[num_dims * sizeof(uint32_t) + 1];
        // std::cout << "number of dimension = " << +num_dims << ", number of levels = " << num_levels << std::endl;
    }
    // exit(0);
    
    std::vector<Type *> reconstructed_vars(n_vars, NULL);

    std::vector<double> var_eb_vec(n_vars, 1);
    for(int i=0; i<n_vars; i++){
        auto reconstructed_data = reconstructors[i].progressive_reconstruct(var_eb_vec[i], -1);
        MGARD::print_statistics(vars_vec[i].data(), reconstructed_data, num_elements);
        MGARD::writefile("reconstructed.dat", reconstructed_data, num_elements);
        reconstructed_vars[i] = reconstructed_data;
    }

}

}
#endif