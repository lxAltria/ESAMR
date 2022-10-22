#ifndef _MDR_ADAPTIVE_RECONSTRUCTOR_HPP
#define _MDR_ADAPTIVE_RECONSTRUCTOR_HPP

#include "ReconstructorInterface.hpp"
#include "Decomposer/Decomposer.hpp"
#include "Interleaver/Interleaver.hpp"
#include "BitplaneEncoder/BitplaneEncoder.hpp"
#include "Retriever/Retriever.hpp"
#include "ErrorEstimator/ErrorEstimator.hpp"
#include "ErrorCollector/ErrorCollector.hpp"
#include "SizeInterpreter/SizeInterpreter.hpp"
#include "LosslessCompressor/LevelCompressor.hpp"
#include "RefactorUtils.hpp"
#include <unordered_map>

namespace MDR {
    // a decomposition-based scientific data reconstructor: inverse operator of composed refactor
    template<class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class SizeInterpreter, class ErrorEstimator, class Retriever>
    class AdaptiveReconstructor : public concepts::ReconstructorInterface<T> {
    public:
        AdaptiveReconstructor(Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, SizeInterpreter interpreter, Retriever retriever)
            : decomposer(decomposer), interleaver(interleaver), encoder(encoder), compressor(compressor), interpreter(interpreter), retriever(retriever){
                load_metadata();
                uint32_t num_elements = 1;
                for(const auto& dim:dims){
                    num_elements *= dim;
                }
                for(int l=0; l<num_levels; l++){
                    level_block_num_bitplanes.push_back(std::vector<uint8_t>(num_blocks[0] * num_blocks[1] * num_blocks[2], 0)); 
                }
                for(int l=0; l<num_levels; l++){
                    const auto& num_agg_blocks = level_num_agg_blocks[l];
                    level_block_num_segments.push_back(std::vector<uint8_t>(num_agg_blocks[0] * num_agg_blocks[1] * num_agg_blocks[2], 0)); 
                }
                converged = std::vector<bool>(num_blocks[0] * num_blocks[1] * num_blocks[2], false);
                data = std::vector<T>(num_elements, 0);
            }

        // reconstruct data from encoded streams
        T * reconstruct(double tolerance){
            // use squared error for convergence tests
            auto prev_level_block_num_segments(level_block_num_segments);
            auto prev_level_block_num_bitplanes(level_block_num_bitplanes);
            // test: retrieve all data
            // TODO: interpret retrieve sizes for each block
            // return retrieved_blocks: retrieved_blocks[l] means required block ids in level l
            //        & level_block_num_bitplanes
            int target_level = 4;
            printf("interpret retrieved_blocks\n");
            std::vector<std::vector<int>> level_retrieved_blocks;
            {
                for(int l=0; l<=target_level; l++){
                    std::vector<int> retrieved_blocks;
                    auto aggregation_granularity = level_aggregation_granularity[l];
                    // compute number of aggregations based on aggregation granularity
                    // agg_nx, agg_ny, agg_nz: number of aggregation blocks along each dimension
                    int agg_nx = (num_blocks[0] < aggregation_granularity) ? num_blocks[0] : aggregation_granularity;
                    int agg_ny = (num_blocks[1] < aggregation_granularity) ? num_blocks[1] : aggregation_granularity;
                    int agg_nz = (num_blocks[2] < aggregation_granularity) ? num_blocks[2] : aggregation_granularity;
                    for(int i=0; i<agg_nx*agg_ny*agg_nz; i++){
                        retrieved_blocks.push_back(i);
                        level_block_num_segments[l][i] = level_block_merge_counts[l][i].size();
                    }
                    for(int i=0; i<num_blocks[0] * num_blocks[1] * num_blocks[2]; i++){
                        level_block_num_bitplanes[l][i] = 32;                        
                    }
                    level_retrieved_blocks.push_back(retrieved_blocks);
                }
            }
            printf("prepare level components\n");
            // prepare level components
            auto total_num_blocks = num_blocks[0] * num_blocks[1] * num_blocks[2];
            level_block_components.clear();
            std::vector<uint8_t *> decompressed_segments;
            for(int l=0; l<=target_level; l++){
                auto block_components = std::vector<std::vector<const uint8_t*>>();
                for(int b=0; b<total_num_blocks; b++){
                    block_components.push_back(std::vector<const uint8_t*>());
                }
                const auto& aggregation_block_map = level_aggregation_block_map[l];
                const auto& agg_block_sizes = level_agg_block_sizes[l];
                const auto& agg_block_bp_sizes = level_agg_block_bp_sizes[l];
                const auto& block_merge_counts = level_block_merge_counts[l];
                const auto& retrieved_blocks = level_retrieved_blocks[l];
                for(int b=0; b<retrieved_blocks.size(); b++){
                    // retrieve one aggregated block of data with multiple precision segment
                    int agg_block_id = retrieved_blocks[b];
                    uint8_t prev_segments = prev_level_block_num_segments[l][agg_block_id];
                    uint8_t curr_segments = level_block_num_segments[l][agg_block_id];
                    // uint8_t num_segments = curr_segments - prev_segments;
                    auto aggregated_level_segments = retriever.retrieve_level_segments(l, agg_block_id, level_agg_block_sizes[l][agg_block_id], prev_segments, curr_segments); 
                    // redistributed to the corresponding blocks
                    int bp_offset = 0;
                    for(int i=0; i<prev_segments; i++){
                        bp_offset += block_merge_counts[agg_block_id][i];
                    }
                    std::cout << "prev_segments = " << +prev_segments << ", curr_segments = " << +curr_segments << std::endl; 
                    for(int i=prev_segments; i<curr_segments; i++){
                        // get id for blocks belonging to this aggregation block
                        const auto& global_blocks = aggregation_block_map.at(agg_block_id);
                        const uint8_t * compressed_data_pos = aggregated_level_segments[i];
                        // deal with each bitplane (which spans all blocks in the aggregation block)
                        // format: b1p1 b2p1 ... bnp1, b1p2 b2p2 ... bnp2, ...
                        for(int bp=0; bp<block_merge_counts[agg_block_id][i]; bp++){
                            uint8_t * precision_segment = NULL;
                            if(l == 0){
                                printf("level = %d, agg_block = %d, segment = %d, bitplane = %d\n", l, agg_block_id, i, bp + bp_offset);
                                printf("agg_block_bp_size = %d\n", agg_block_bp_sizes[agg_block_id][bp + bp_offset]);
                            }
                            // printf("offset = %ld\n", compressed_data_pos - aggregated_level_segments[i]);
                            // for(int i=0; i<agg_block_bp_sizes[agg_block_id][bp + bp_offset]; i++){
                            //     std::cout << +compressed_data_pos[i] << " ";
                            // }
                            // std::cout << std::endl;
                            auto decompressed_size = ZSTD::decompress(compressed_data_pos, agg_block_bp_sizes[agg_block_id][bp + bp_offset], &precision_segment);
                            // printf("decompressed_size = %d\n", decompressed_size);
                            // for(int i=0; i<decompressed_size; i++){
                            //     std::cout << +precision_segment[i] << " ";
                            // }
                            // std::cout << std::endl;
                            decompressed_segments.push_back(precision_segment);
                            // prepare each precision fragment
                            const uint8_t * bitplane_pos = precision_segment;
                            for(int j=0; j<global_blocks.size(); j++){
                                auto block_id = global_blocks[j];
                                auto block_num_element = block_level_elements[block_id][l];  
                                // TODO: change to adaptive computation with respect to T_stream
                                auto block_stream_size = ((block_num_element - 1) / (UINT8_BITS*4) + 1) * sizeof(uint32_t);                 
                                // printf("block_id = %d, block_num_element = %d, block_stream_size = %d\n", block_id, block_num_element, block_stream_size);
                                if(!converged[block_id]){
                                    block_components[block_id].push_back(bitplane_pos);
                                }
                                bitplane_pos += block_stream_size;
                            }
                            compressed_data_pos += agg_block_bp_sizes[agg_block_id][bp + bp_offset];
                        }
                        bp_offset += block_merge_counts[agg_block_id][i];
                    }
                    // exit(0);                    
                }
                level_block_components.push_back(block_components);
            }
            bool success = reconstruct(target_level, prev_level_block_num_bitplanes);
            printf("free decompressed data\n");
            for(int i=0; i<decompressed_segments.size(); i++){
                free(decompressed_segments[i]);
            }
            return data.data();
        }

        // reconstruct progressively based on available data
        // Note: the returned data may not be used directly
        T * progressive_reconstruct(double tolerance){
            std::vector<T> cur_data(data);
            reconstruct(tolerance);
            for(int i=0; i<data.size(); i++){
                data[i] += cur_data[i];
            }                
            return data.data();
        }

        void load_metadata(){
            uint8_t * metadata = retriever.load_metadata();
            uint8_t const * metadata_pos = metadata;
            uint8_t num_dims = *(metadata_pos ++);
            deserialize(metadata_pos, num_dims, dims);
            this->num_bitplanes = *(metadata_pos ++);
            this->block_size = *(metadata_pos ++);
            num_levels = *(metadata_pos ++);
            std::vector<std::vector<int>> level_merge_counts;
            std::vector<std::vector<std::vector<double>>> level_squared_errors;
            std::vector<std::vector<uint32_t>> level_sizes;
            deserialize(metadata_pos, num_levels, level_max_exp);
            deserialize(metadata_pos, num_levels, level_merge_counts);
            deserialize(metadata_pos, num_levels, level_squared_errors);
            // deserialize(metadata_pos, num_levels, level_sizes);
            deserialize(metadata_pos, num_levels, level_num);
            deserialize(metadata_pos, num_levels, level_agg_block_bp_sizes);
            // read done
            num_blocks = std::vector<uint32_t>(dims.size());
            for(int i=0; i<dims.size(); i++){
                num_blocks[i] = (dims[i] - 1) / block_size + 1;
            }
            strides = std::vector<uint32_t>(dims.size());
            uint32_t stride = 1;
            for(int i=dims.size()-1; i>=0; i--){
                strides[i] = stride;
                stride *= dims[i];
            } 
            auto total_num_blocks = num_blocks[0] * num_blocks[1] * num_blocks[2];
            // print_vec("merge counts", level_merge_counts);
            // print_vec("sizes", level_sizes);
            // print_vec("merged errors", level_squared_errors);
            // compute level_aggregation_granularity and level_block_merge_counts
            level_aggregation_granularity.clear();
            level_block_merge_counts.clear();
            for(int i=0; i<num_levels; i++){
                int bp_num = 0;
                int agg_block_id = 0;
                int block_sq_err_id = 0;
                std::vector<std::vector<int>> block_merge_counts;
                std::vector<int> merge_counts;
                for(int j=0; j<level_merge_counts[i].size(); j++){
                    bp_num += level_merge_counts[i][j];
                    merge_counts.push_back(level_merge_counts[i][j]);
                    if(bp_num % num_bitplanes == 0){
                        // end of one block
                        block_merge_counts.push_back(merge_counts);
                        merge_counts.clear();
                        agg_block_id ++;
                    }
                }
                level_aggregation_granularity.push_back(sqrt(agg_block_id));
                level_block_merge_counts.push_back(block_merge_counts);
            }
            build_block_aggregation_map();
            // compute level_num_agg_blocks
            for(int l=0; l<num_levels; l++){
                auto aggregation_granularity = level_aggregation_granularity[l];
                int agg_nx = (num_blocks[0] < aggregation_granularity) ? num_blocks[0] : aggregation_granularity;
                int agg_ny = (num_blocks[1] < aggregation_granularity) ? num_blocks[1] : aggregation_granularity;
                int agg_nz = (num_blocks[2] < aggregation_granularity) ? num_blocks[2] : aggregation_granularity;
                std::vector<uint32_t> num_agg_blocks;
                num_agg_blocks.push_back(agg_nx);
                num_agg_blocks.push_back(agg_ny);
                num_agg_blocks.push_back(agg_nz);
                level_num_agg_blocks.push_back(num_agg_blocks);
            }
            // compute level_block_squared_errors
            level_block_squared_errors.clear();
            for(int l=0; l<num_levels; l++){
                std::vector<std::vector<double>> block_squared_errors;
                for(int b=0; b<total_num_blocks; b++){
                    block_squared_errors.push_back(std::vector<double>());
                }
                level_block_squared_errors.push_back(block_squared_errors);
            }
            for(int l=0; l<num_levels; l++){
                const auto& squared_errors = level_squared_errors[l];
                int precision_index = 0;
                for(int i=0; i<level_block_merge_counts[l].size(); i++){
                    // the i-th agg_block
                    const auto& block_ids = level_aggregation_block_map[l][i];
                    for(int j=0; j<block_ids.size(); j++){
                        auto id = block_ids[j];
                        // the 0-th precion
                        level_block_squared_errors[l][id].push_back(squared_errors[precision_index][j]);
                        for(int k=0; k<level_block_merge_counts[l][i].size(); k++){
                            // the k-th precision
                            level_block_squared_errors[l][id].push_back(squared_errors[precision_index + 1 + k][j]);
                        }
                    }
                    precision_index += level_block_merge_counts[l][i].size() + 1;
                }
                
            }
            // compute level_agg_block_sizes
            level_agg_block_sizes.clear();
            for(int l=0; l<num_levels; l++){
                std::vector<std::vector<uint32_t>> agg_block_sizes;
                for(int b=0; b<level_block_merge_counts[l].size(); b++){
                    const auto& merge_counts = level_block_merge_counts[l][b];
                    const auto& bp_sizes = level_agg_block_bp_sizes[l][b];
                    std::vector<uint32_t> sizes;
                    int index = 0;
                    for(int i=0; i<merge_counts.size(); i++){
                        size_t segment_size = 0;
                        for(int j=0; j<merge_counts[i]; j++){
                            segment_size += bp_sizes[index];
                            index ++;
                        }
                        sizes.push_back(segment_size);
                    }
                    agg_block_sizes.push_back(sizes);
                }
                level_agg_block_sizes.push_back(agg_block_sizes);
            }
            // print_vec(level_aggregation_granularity);
            // for(int l=0; l<num_levels; l++){
            //     print_vec("level_agg_block_sizes", level_agg_block_sizes[l]);
            //     print_vec("level_agg_block_bp_sizes", level_agg_block_bp_sizes[l]);
            //     print_vec("level_block_merge_counts", level_block_merge_counts[l]);                
            //     print_vec("level_block_squared_errors", level_block_squared_errors[l]);                
            // }
            // compute block_level_dims and block_level_elements
            auto nx = num_blocks[0];
            auto ny = num_blocks[1];
            auto nz = num_blocks[2];
            std::vector<uint32_t> block_dims(3, block_size);
            for(int i=0; i<nx; i++){
                block_dims[0] = (i < nx - 1) ? block_size : (dims[0] - i*block_size);
                for(int j=0; j<ny; j++){
                    block_dims[1] = (j < ny - 1) ? block_size : (dims[1] - j*block_size);
                    for(int k=0; k<nz; k++){
                        block_dims[2] = (k < nz - 1) ? block_size : (dims[2] - k*block_size);
                        auto level_dims = compute_level_dims(block_dims, num_levels - 1);
                        auto level_elements = compute_level_elements(level_dims, num_levels - 1);
                        block_level_dims.push_back(level_dims);
                        block_level_elements.push_back(level_elements);
                    }
                }
            }
            free(metadata);
        }

        const std::vector<uint32_t>& get_dimensions(){
            return dims;
        }

        ~AdaptiveReconstructor(){}

        void print() const {
            std::cout << "Adaptive reconstructor with the following components." << std::endl;
            std::cout << "Decomposer: "; decomposer.print();
            std::cout << "Interleaver: "; interleaver.print();
            std::cout << "Encoder: "; encoder.print();
            std::cout << "SizeInterpreter: "; interpreter.print();
            std::cout << "Retriever: "; retriever.print();
        }
    private:
        struct AggregationBlockInfo{
            int agg_block_id;       // id of aggregation block
            int agg_block_offset;   // offset inside the aggregation block
            AggregationBlockInfo(){}
            AggregationBlockInfo(int id, int offset){
                agg_block_id = id;
                agg_block_offset = offset;
            }
        };

        void build_block_aggregation_map(){
            level_block_aggregation_map.clear();
            for(int l=0; l<level_max_exp.size(); l++){
                std::unordered_map<int, AggregationBlockInfo> block_aggregation_map;
                std::unordered_map<int, std::vector<int>> aggregation_block_map;
                auto aggregation_granularity = level_aggregation_granularity[l];
                // compute number of aggregations based on aggregation granularity
                // agg_nx, agg_ny, agg_nz: number of aggregation blocks along each dimension
                int agg_nx = (num_blocks[0] < aggregation_granularity) ? num_blocks[0] : aggregation_granularity;
                int agg_ny = (num_blocks[1] < aggregation_granularity) ? num_blocks[1] : aggregation_granularity;
                int agg_nz = (num_blocks[2] < aggregation_granularity) ? num_blocks[2] : aggregation_granularity;
                // agg_block_size_x, agg_block_size_y, agg_block_size_z: number of data blocks in each aggregation block
                int agg_block_size_x = num_blocks[0] / agg_nx;
                int agg_block_size_y = num_blocks[1] / agg_ny;
                int agg_block_size_z = num_blocks[2] / agg_nz;
                // residual for computing block size
                int residual_x = num_blocks[0] % agg_nx;
                int residual_y = num_blocks[1] % agg_ny;
                int residual_z = num_blocks[2] % agg_nz;
                // iterate each aggregation area
                for(int i=0; i<agg_nx; i++){
                    int agg_actual_block_size_nx = (i < residual_x) ? (agg_block_size_x + 1) : agg_block_size_x;
                    int block_offset_x = (i < residual_x) ? (i * (agg_block_size_x + 1)) : (i*agg_block_size_x + residual_x);
                    for(int j=0; j<agg_ny; j++){
                        int agg_actual_block_size_ny = (j < residual_y) ? (agg_block_size_y + 1) : agg_block_size_y;
                        int block_offset_y = (j < residual_y) ? (j * (agg_block_size_y + 1)) : (j*agg_block_size_y + residual_y);
                        for(int k=0; k<agg_nz; k++){
                            int agg_actual_block_size_nz = (k < residual_z) ? (agg_block_size_z + 1) : agg_block_size_z;
                            int block_offset_z = (k < residual_z) ? (k * (agg_block_size_z + 1)) : (k*agg_block_size_z + residual_z);
                            int num_aggregation = agg_actual_block_size_nx * agg_actual_block_size_ny * agg_actual_block_size_nz;
                            // iterate data in one aggregation area which has aggregation_size^3 blocks
                            // buffer is used to store the encoded data
                            int aggregation_block_id = i*agg_ny*agg_nz + j*agg_nz + k;
                            int aggregation_id = 0;
                            std::vector<int> aggregation_blocks;
                            for(int ii=0; ii<agg_actual_block_size_nx; ii++){
                                for(int jj=0; jj<agg_actual_block_size_ny; jj++){
                                    for(int kk=0; kk<agg_actual_block_size_nz; kk++){
                                        // deal with data in one block 
                                        int block_id = (block_offset_x + ii) * num_blocks[1] * num_blocks[2] + (block_offset_y + jj) * num_blocks[2] + (block_offset_z + kk);
                                        block_aggregation_map.insert({block_id, AggregationBlockInfo(aggregation_block_id, aggregation_id)});
                                        aggregation_blocks.push_back(block_id);
                                        aggregation_id ++;
                                    }
                                }
                            }
                            aggregation_block_map.insert({aggregation_block_id, aggregation_blocks});
                        }
                    }
                }
                level_block_aggregation_map.push_back(block_aggregation_map);
                level_aggregation_block_map.push_back(aggregation_block_map);
            }
            // print map
            // for(int l=0; l<level_block_aggregation_map.size(); l++){
            //     const auto& map = level_block_aggregation_map[l];
            //     for(const auto& iter:map){
            //         std::cout << iter.first << " (" << iter.second.agg_block_id << ", " << iter.second.agg_block_offset << ")\n";
            //     }
            // }
        }

        template<class T1>
        std::vector<std::vector<T1>> extract_block_info(const std::vector<std::vector<std::vector<T1>>>& level_block_info, int block_id){
            std::vector<std::vector<T1>> block_info;
            for(int i=0; i<level_block_info.size(); i++){
                block_info.push_back(level_block_info[i][block_id]);
            }
            return block_info;
        }

        // reconstruct data based on given level components
        bool reconstruct(uint8_t target_level, const std::vector<std::vector<uint8_t>>& prev_level_block_num_bitplanes, bool progressive=true){
            memset(data.data(), 0, data.size() * sizeof(T));
            std::cout << "target_level = " << +target_level << ", block_size = " << block_size << std::endl;
            print_vec(num_blocks);
            auto nx = num_blocks[0];
            auto ny = num_blocks[1];
            auto nz = num_blocks[2];
            uint32_t total_num_blocks = nx * ny * nz;
            std::vector<uint32_t> block_dims(3, block_size);
            uint32_t max_num_block_elements = block_size * block_size * block_size;
            auto max_level_dims = compute_level_dims(block_dims, target_level);
            auto max_level_elements = compute_level_elements(max_level_dims, target_level);
            std::vector<std::vector<uint32_t>> level_block_elements;
            // decompose and interleave
            // the buffer store data in a blockwise then levelwise fashtion
            // i.e.: |B1L1, B1L2, ... B2L1, B2L2 ...|
            // stride between different block is max_num_block_elements
            T * buffer = (T *) malloc(max_num_block_elements * total_num_blocks * sizeof(T));
            std::cout << "max_num_block_elements = " << max_num_block_elements << ", total_num_blocks = " << total_num_blocks << std::endl;
            int block_id = 0;
            std::vector<uint32_t> dims_dummy(dims.size(), 0);
            T * buffer_pos = buffer;
            T * x_data_pos = data.data();
            for(int i=0; i<nx; i++){
                block_dims[0] = (i < nx - 1) ? block_size : (dims[0] - i*block_size);
                T * y_data_pos = x_data_pos;
                for(int j=0; j<ny; j++){
                    block_dims[1] = (j < ny - 1) ? block_size : (dims[1] - j*block_size);
                    T * z_data_pos = y_data_pos;
                    for(int k=0; k<nz; k++){
                        block_dims[2] = (k < nz - 1) ? block_size : (dims[2] - k*block_size);
                        if(!converged[block_id]){
                            // decompose and interleave data in each block
                            // print_vec(block_dims);
                            const auto& level_dims = block_level_dims[block_id];
                            const auto& level_elements = block_level_elements[block_id];
                            const auto& reconstruct_dimensions = level_dims[target_level];
                            for(int l=0; l<=target_level; l++){
                                const auto& components = level_block_components[l][block_id];
                                if(components.size() > 0){
                                    auto prev_num_bp = prev_level_block_num_bitplanes[l][block_id];
                                    auto curr_num_bp = level_block_num_bitplanes[l][block_id];
                                    int agg_block_id = level_block_aggregation_map[l][block_id].agg_block_id;
                                    // debug
                                    // if(block_id == 1){
                                    //     std::cout << "l = " << l << ", block_id = " << block_id << ", agg_block_id = " << agg_block_id << std::endl;
                                    //     std::cout << "level_elements = " << level_elements[l] << ", prev_num_bp = " << +prev_num_bp << ", curr_num_bp = " << +curr_num_bp << std::endl;
                                    // }
                                    T * level_decoded_data = encoder.progressive_decode(components, level_elements[l], level_max_exp[l], prev_num_bp, curr_num_bp - prev_num_bp, l);
                                    const std::vector<uint32_t>& prev_dims = (l == 0) ? dims_dummy : level_dims[l - 1];
                                    interleaver.reposition(level_decoded_data, reconstruct_dimensions, level_dims[l], prev_dims, z_data_pos, this->strides);
                                    free(level_decoded_data);
                                }
                            }
                            printf("recompose block %d\n", block_id);
                            // std::cout << reconstruct_dimensions[0] << " " << reconstruct_dimensions[1] << " " << reconstruct_dimensions[2] << std::endl;
                            // printf("target_level = %d\n", target_level);
                            decomposer.recompose(z_data_pos, reconstruct_dimensions, target_level, this->strides);
                        }
                        block_id ++;
                        z_data_pos += block_size;
                    }
                    y_data_pos += block_size * dims[2];
                }
                x_data_pos += block_size * dims[1] * dims[2];
            }
            return true;
        }

        Decomposer decomposer;
        Interleaver interleaver;
        Encoder encoder;
        SizeInterpreter interpreter;
        Retriever retriever;
        Compressor compressor;
        std::vector<T> data;
        std::vector<uint32_t> dims;
        std::vector<int> level_aggregation_granularity;
        std::vector<int> level_max_exp;
        std::vector<std::vector<std::vector<int>>> level_block_merge_counts;
        std::vector<std::vector<std::vector<double>>> level_block_squared_errors;
        std::vector<std::vector<std::vector<const uint8_t*>>> level_block_components;
        std::vector<std::vector<std::vector<uint32_t>>> level_agg_block_sizes;
        std::vector<std::vector<std::vector<uint32_t>>> level_agg_block_bp_sizes;
        std::vector<bool> converged;
        std::vector<uint32_t> level_num;
        std::vector<uint32_t> num_blocks;
        std::vector<std::vector<uint32_t>> level_num_agg_blocks;
        std::vector<std::vector<uint8_t>> level_block_num_bitplanes; // L x B
        std::vector<std::vector<uint8_t>> level_block_num_segments; // sum(aggB)
        std::vector<std::vector<std::vector<uint32_t>>> block_level_dims;
        std::vector<std::vector<uint32_t>> block_level_elements;
        std::vector<uint32_t> strides;
        uint8_t num_levels = 0;
        int block_size = 0;
        int num_bitplanes = 0;
        int num_elements = 0;
        int current_level = 0;
        std::vector<std::unordered_map<int, AggregationBlockInfo>> level_block_aggregation_map; // map global block to aggregation block
        std::vector<std::unordered_map<int, std::vector<int>>> level_aggregation_block_map;     // map aggregation block to global blocks
    };
}
#endif

