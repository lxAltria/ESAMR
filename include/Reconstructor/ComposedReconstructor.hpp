#ifndef _MDR_COMPOSED_RECONSTRUCTOR_HPP
#define _MDR_COMPOSED_RECONSTRUCTOR_HPP

#include "ReconstructorInterface.hpp"
#include "Decomposer/Decomposer.hpp"
#include "Interleaver/Interleaver.hpp"
#include "BitplaneEncoder/BitplaneEncoder.hpp"
#include "Retriever/Retriever.hpp"
#include "ErrorEstimator/ErrorEstimator.hpp"
#include "LosslessCompressor/LosslessCompressor.hpp"
#include "RefactorUtils.hpp"

namespace MDR {
    // a decomposition-based scientific data reconstructor: inverse operator of composed refactor
    template<class T, class Decomposer, class Interleaver, class Encoder, class Retriever>
    class ComposedReconstructor : public concepts::ReconstructorInterface<T> {
    public:
        ComposedReconstructor(Decomposer decomposer, Interleaver interleaver, Encoder encoder, Retriever retriever)
            : decomposer(decomposer), interleaver(interleaver), encoder(encoder), retriever(retriever){}

        T * reconstruct(uint8_t const * refactored_data, double tolerance){
            level_num_bitplanes.clear();
            uint32_t retrieve_size = retriever.interpret_size(level_sizes, level_errors, level_order, tolerance, level_num_bitplanes);
            load_level_components(refactored_data, retrieve_size);
            // check whether to reconstruct to full resolution
            uint8_t target_level = level_error_bounds.size() - 1;
            int skipped_level = 0;
            for(int i=0; i<=target_level; i++){
                if(level_num_bitplanes[target_level - i] != 0){
                    skipped_level = i;
                    break;
                }
            }
            target_level -= skipped_level;

            if(reconstruct(target_level)) return data.data();
            else{
                std::cerr << "Recontruct unsuccessful, return NULL pointer" << std::endl;
                return NULL;
            }
        }

        void load_metadata(uint8_t const * metadata){
            uint8_t const * metadata_pos = metadata;
            uint8_t num_dims = *(metadata_pos ++);
            deserialize(metadata_pos, num_dims, dimensions);
            uint8_t num_levels = *(metadata_pos ++);
            deserialize(metadata_pos, num_levels, level_error_bounds);
            deserialize(metadata_pos, num_levels, level_errors);
            deserialize(metadata_pos, num_levels, level_sizes);
            uint32_t num_level_orders = *reinterpret_cast<const uint32_t*>(metadata_pos);
            metadata_pos += sizeof(uint32_t);
            deserialize(metadata_pos, num_level_orders, level_order);
        }

        ~ComposedReconstructor(){}

        void print() const {
            std::cout << "Composed reconstructor with the following components." << std::endl;
            std::cout << "Decomposer: "; decomposer.print();
            std::cout << "Interleaver: "; interleaver.print();
            std::cout << "Encoder: "; encoder.print();
            std::cout << "Retriever: "; retriever.print();
        }
    private:

        bool load_level_components(uint8_t const * refactored_data, uint32_t retrieve_size){
            level_components.clear();
            for(int i=0; i<level_sizes.size(); i++){
                level_components.push_back(std::vector<const uint8_t*>());
            }
            const uint8_t * refactored_data_pos = refactored_data;
            std::vector<int> index(level_sizes.size(), 0);
            int count = 0;
            while(refactored_data_pos - refactored_data < retrieve_size){
                int level = level_order[count ++];
                int bitplane_index = index[level];
                level_components[level].push_back(refactored_data_pos);
                refactored_data_pos += level_sizes[level][bitplane_index];
                index[level] ++;
            }
            return true;
        }

        bool reconstruct(uint8_t target_level){
            auto level_dims = compute_level_dims(dimensions, target_level);
            auto reconstruct_dimensions = level_dims[target_level];
            uint32_t num_elements = 1;
            for(const auto& dim:reconstruct_dimensions){
                num_elements *= dim;
            }
            data.clear();
            data = std::vector<T>(num_elements, 0);

            auto lossless_compressor = ZSTD();
            auto level_elements = compute_level_elements(level_dims, target_level);
            std::vector<uint32_t> dims_dummy(reconstruct_dimensions.size(), 0);
            for(int i=0; i<=target_level; i++){
                std::vector<const uint8_t*> decompressed_level_compenents;
                std::vector<uint8_t*> lossless_decompressed;
                for(int j=0; j<level_num_bitplanes[i]; j++){
                    uint8_t * decompressed = NULL;
                    auto decompressed_size = lossless_compressor.decompress(level_components[i][j], level_sizes[i][j], &decompressed);
                    decompressed_level_compenents.push_back(decompressed);
                    if(decompressed != level_components[i][j]){
                        lossless_decompressed.push_back(decompressed);                    
                    }
                }
                int level_exp = 0;
                frexp(level_error_bounds[i], &level_exp);
                auto level_decoded_data = encoder.decode(decompressed_level_compenents, level_elements[i], level_exp, level_num_bitplanes[i]);
                for(int j=0; j<lossless_decompressed.size(); j++){
                    free(lossless_decompressed[j]);
                }

                const std::vector<uint32_t>& prev_dims = (i == 0) ? dims_dummy : level_dims[i - 1];
                interleaver.reposition(level_decoded_data, reconstruct_dimensions, level_dims[i], prev_dims, data.data());
                free(level_decoded_data);
            }
            decomposer.recompose(data.data(), reconstruct_dimensions, target_level);
            return true;
        }

        Decomposer decomposer;
        Interleaver interleaver;
        Encoder encoder;
        Retriever retriever;
        std::vector<T> data;
        std::vector<uint32_t> dimensions;
        std::vector<T> level_error_bounds;
        std::vector<uint8_t> level_num_bitplanes;
        std::vector<std::vector<const uint8_t*>> level_components;
        std::vector<std::vector<uint32_t>> level_sizes;
        std::vector<std::vector<double>> level_errors;
        std::vector<uint8_t> level_order;
    };
}
#endif

