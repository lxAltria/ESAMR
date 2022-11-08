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
#include <set>
#include <queue>

namespace MDR {
    // a decomposition-based scientific data reconstructor: inverse operator of composed refactor
    template<class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class SizeInterpreter, class ErrorEstimator, class Retriever>
    class AdaptiveReconstructor : public concepts::ReconstructorInterface<T> {
    public:
        AdaptiveReconstructor(Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, SizeInterpreter interpreter, Retriever retriever)
            : decomposer(decomposer), interleaver(interleaver), encoder(encoder), compressor(compressor), interpreter(interpreter), retriever(retriever){
                load_metadata();
                num_elements = 1;
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
                // TODO: use estimator as input parameter
                error_estimator = SNormErrorEstimator<T>(dims.size(), level_max_exp.size() - 1, 0.0);
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
            auto total_num_blocks = num_blocks[0] * num_blocks[1] * num_blocks[2];
            int target_level = level_max_exp.size() - 1;
            printf("interpret retrieved_blocks\n");
            size_t retrieved_size = 0;
            std::vector<std::vector<int>> level_retrieved_blocks;
            {
                // using a variational greedy based method
                // the level-block components build a tree structure
                // using greedy algorithm for each block
                // compare adjustacent level by combining the impact of all related blocks
                // for(int l=0; l<=target_level; l++){
                //     print_vec("level_block_squared_errors", level_block_squared_errors[l]);
                // }
                // for(int l=0; l<=target_level; l++){
                //     print_vec("level_agg_block_sizes", level_agg_block_sizes[l]);
                // }
                // build the hierarchical tree for aggregation blocks
                // using set for construction and then convert to vectors
                std::vector<std::vector<std::set<int>>> agg_block_hierarchy_set;
                for(int l=0; l<target_level; l++){
                    std::vector<std::set<int>> hierarchy;
                    for(int i=0; i<level_aggregation_block_map[l].size(); i++){
                        hierarchy.push_back(std::set<int>());
                    }
                    agg_block_hierarchy_set.push_back(hierarchy);
                }
                for(int i=0; i<total_num_blocks; i++){
                    for(int l=0; l<target_level; l++){
                        auto curr_agg_block_id = level_block_aggregation_map[l][i].agg_block_id;
                        auto child_agg_block_id = level_block_aggregation_map[l+1][i].agg_block_id;
                        agg_block_hierarchy_set[l][curr_agg_block_id].insert(child_agg_block_id);
                    }
                }
                // agg_block_hierarchy: [l][agg_block_id] -> agg_block_id in l+1
                std::vector<std::vector<std::vector<int>>> agg_block_hierarchy;
                for(int l=0; l<target_level; l++){
                    std::vector<std::vector<int>> hierarchy;
                    for(int i=0; i<agg_block_hierarchy_set[l].size(); i++){
                        hierarchy.push_back(std::vector<int>(agg_block_hierarchy_set[l][i].begin(), agg_block_hierarchy_set[l][i].end()));
                    }
                    agg_block_hierarchy.push_back(hierarchy);
                }
                // reversed_agg_block_hierarchy: [l][agg_block_id] -> agg_block_id in l-1
                std::vector<std::vector<int>> reversed_agg_block_hierarchy;
                // place holder for level 0
                reversed_agg_block_hierarchy.push_back(std::vector<int>());
                for(int l=1; l<=target_level; l++){
                    reversed_agg_block_hierarchy.push_back(std::vector<int>(level_aggregation_block_map[l].size()));
                }
                for(int l=0; l<target_level; l++){
                    for(int i=0; i<agg_block_hierarchy[l].size(); i++){
                        for(int j=0; j<agg_block_hierarchy[l][i].size(); j++){
                            int agg_block_id = agg_block_hierarchy[l][i][j];
                            reversed_agg_block_hierarchy[l+1][agg_block_id] = i;
                        }
                    }
                    
                }
                // for(int l=0; l<target_level; l++){
                //     printf("agg_block_hierarchy in level %d\n", l);
                //     for(int i=0; i<agg_block_hierarchy[l].size(); i++){
                //         for(const auto& val:agg_block_hierarchy[l][i]){
                //             std::cout << val << " ";
                //         }
                //         std::cout << std::endl;
                //     }
                // }
                // print_vec("reversed_agg_block_hierarchy", reversed_agg_block_hierarchy);
                // compute retrieved block using a greedy algorithm
                std::vector<std::vector<double>> level_block_error_gain;
                std::vector<std::vector<uint32_t>> level_block_fetch_size;
                std::vector<std::vector<bool>> level_block_fetch_option; // true for lower level, false for current level
                for(int l=0; l<=target_level; l++){
                    int level_num_agg_blocks = level_aggregation_block_map[l].size();
                    level_block_error_gain.push_back(std::vector<double>(level_num_agg_blocks, 0));
                    level_block_fetch_size.push_back(std::vector<uint32_t>(level_num_agg_blocks, 0));
                    level_block_fetch_option.push_back(std::vector<bool>(level_num_agg_blocks, false));
                    if(l == target_level){
                        // init statistics for the last level
                        for(int i=0; i<level_aggregation_block_map[l].size(); i++){
                            auto curr_num_segments = level_block_num_segments[l][i];
                            level_block_fetch_size[l][i] = level_agg_block_sizes[l][i][curr_num_segments];
                            // iterate through all blocks in the aggregation
                            for(const auto& block_id : level_aggregation_block_map[l][i]){
                                if(!converged[block_id]){
                                    level_block_error_gain[l][i] += error_estimator.estimate_error_gain(0, level_block_squared_errors[l][block_id][curr_num_segments], level_block_squared_errors[l][block_id][curr_num_segments+1], l);
                                }
                            }
                            level_block_fetch_option[l][i] = false;
                        }
                    }
                }
                for(int l=target_level-1; l>=0; l--){
                    // iterate through all aggregation block in the level
                    for(int i=0; i<level_aggregation_block_map[l].size(); i++){
                        // check the efficiency of segment in current level
                        double error_gain_curr_level = 0;
                        uint32_t size_curr_level = 0;
                        double efficiency_curr_level = compute_efficiency(l, i, level_block_num_segments[l][i], level_aggregation_block_map[l][i], error_gain_curr_level, size_curr_level);
                        // check the efficiency of segment in the lower level
                        double error_gain_lower_level = 0;
                        uint32_t size_lower_level = 0;
                        for(int j=0; j<agg_block_hierarchy[l][i].size(); j++){
                            auto block_id = agg_block_hierarchy[l][i][j];
                            error_gain_lower_level += level_block_error_gain[l+1][block_id];
                            size_lower_level += level_block_fetch_size[l+1][block_id];
                        }
                        double efficiency_lower_level = error_gain_lower_level / size_lower_level;
                        if(efficiency_curr_level > efficiency_lower_level){
                            // choose segment in current level
                            level_block_fetch_option[l][i] = false;
                            level_block_fetch_size[l][i] = size_curr_level;
                            level_block_error_gain[l][i] = error_gain_curr_level;
                        }
                        else{
                            // choose lower level
                            level_block_fetch_option[l][i] = true;
                            level_block_fetch_size[l][i] = size_lower_level;
                            level_block_error_gain[l][i] = error_gain_lower_level;
                        }
                    }
                }
                // print_vec("level_block_error_gain", level_block_error_gain);
                // print_vec("level_block_fetch_size", level_block_fetch_size);
                // print_vec("level_block_fetch_option", level_block_fetch_option);
                double estimated_error = 0;
                for(int l=0; l<=target_level; l++){
                    for(int i=0; i<total_num_blocks; i++){
                        auto curr_num_segments = level_block_num_segments[l][i];
                        estimated_error += error_estimator.estimate_error(level_block_squared_errors[l][i][curr_num_segments], l);
                    }                    
                }
                std::vector<std::set<int>> level_retrieved_blocks_set;
                for(int i=0; i<=target_level; i++){
                    level_retrieved_blocks_set.push_back(std::set<int>());
                }
                int count = 0;
                // std::cout << "estimated_error = " << estimated_error << std::endl;
                while(estimated_error > tolerance){
                    if(count > 2){
                        break;
                    }
                    // std::cout << "round " << count << std::endl;
                    count ++;
                    // find the path to fetch options
                    // using a queue of <level, agg_block_id> pair
                    std::queue<std::pair<int, int>> agg_block_queue;
                    agg_block_queue.push({0, 0});
                    while(!agg_block_queue.empty()){
                        // std::cout << "starting processing queue" << std::endl;
                        int level = agg_block_queue.front().first;
                        int agg_block_id = agg_block_queue.front().second;
                        agg_block_queue.pop();
                        if(level_block_fetch_option[level][agg_block_id]){
                            // need to retrieve next level
                            for(int i=0; i<agg_block_hierarchy[level][agg_block_id].size(); i++){
                                agg_block_queue.push({level + 1, agg_block_hierarchy[level][agg_block_id][i]});
                            }
                        }
                        else if(level_block_num_segments[level][agg_block_id] < level_block_merge_counts[level][agg_block_id].size()){
                            // retrieve current level when size is not exceeded
                            auto curr_num_segments = level_block_num_segments[level][agg_block_id];
                            // increment total retrieved size
                            retrieved_size += level_agg_block_sizes[level][agg_block_id][curr_num_segments];
                            // printf("retrieve L%d_B%d_S%d\n", level, agg_block_id, curr_num_segments);
                            level_retrieved_blocks_set[level].insert(agg_block_id);
                            // update level_block_num_segments and level_block_num_bitplanes
                            auto bitplane_count = level_block_merge_counts[level][agg_block_id][curr_num_segments];
                            level_block_num_segments[level][agg_block_id] ++;
                            for(const auto& block_id : level_aggregation_block_map[level][agg_block_id]){
                                level_block_num_bitplanes[level][block_id] += bitplane_count;
                            }
                            // need to update efficiency
                            // printf("update efficiency\n");
                            int l = level;
                            int i = agg_block_id;
                            while(l >= 0){
                                // printf("compute curr_level efficiency\n");
                                // compute efficiency of level l
                                double efficiency_curr_level = 0;
                                double error_gain_curr_level = 0;
                                uint32_t size_curr_level = 0;
                                if(level_block_num_segments[l][i] < level_block_merge_counts[l][i].size()){
                                    efficiency_curr_level = compute_efficiency(l, i, level_block_num_segments[l][i], level_aggregation_block_map[l][i], error_gain_curr_level, size_curr_level);                                    
                                }
                                // printf("compute lower_level efficiency\n");
                                double error_gain_lower_level = 0;
                                uint32_t size_lower_level = 0;
                                // compute efficiency of level l+1
                                if(l < target_level){
                                    for(int j=0; j<agg_block_hierarchy[l][i].size(); j++){
                                        auto agg_block_id = agg_block_hierarchy[l][i][j];
                                        error_gain_lower_level += level_block_error_gain[l+1][agg_block_id];
                                        size_lower_level += level_block_fetch_size[l+1][agg_block_id];
                                    }                                    
                                }
                                if((!error_gain_curr_level) && (!error_gain_lower_level)){
                                    break;
                                }
                                bool select_current_level = true;
                                if(!error_gain_curr_level) select_current_level = false;
                                else if(!error_gain_lower_level) select_current_level = true;
                                else{
                                    double efficiency_lower_level = error_gain_lower_level / size_lower_level;
                                    select_current_level = (efficiency_curr_level > efficiency_lower_level);
                                }
                                if(select_current_level){
                                    // choose segment in current level
                                    level_block_fetch_option[l][i] = false;
                                    level_block_fetch_size[l][i] = size_curr_level;
                                    level_block_error_gain[l][i] = error_gain_curr_level;
                                }
                                else{
                                    // choose lower level
                                    level_block_fetch_option[l][i] = true;
                                    level_block_fetch_size[l][i] = size_lower_level;
                                    level_block_error_gain[l][i] = error_gain_lower_level;
                                }
                                // printf("processing done\n");
                                if(l > 0) i = reversed_agg_block_hierarchy[l][i];
                                l = l - 1;
                                // printf("initializing block %d in level %d\n", i, l);
                            }
                        }
                    }
                }
                for(int l=0; l<=target_level; l++){
                    level_retrieved_blocks.push_back(std::vector<int>(level_retrieved_blocks_set[l].begin(), level_retrieved_blocks_set[l].end()));
                }
            }
            // print_vec("level_retrieved_blocks", level_retrieved_blocks);
            // print_vec("prev_level_block_num_segments", prev_level_block_num_segments);
            // print_vec("level_block_num_segments", level_block_num_segments);
            // print_vec("prev_level_block_num_bitplanes", prev_level_block_num_bitplanes);
            print_vec("level_block_num_bitplanes", level_block_num_bitplanes);
            printf("prepare level components\n");
            // prepare level components
            level_block_components.clear();
            std::vector<uint8_t *> decompressed_segments;
            for(int l=0; l<=target_level; l++){
                auto block_components = std::vector<std::vector<const uint8_t*>>();
                for(int b=0; b<total_num_blocks; b++){
                    block_components.push_back(std::vector<const uint8_t*>());
                }
                if(level_retrieved_blocks[l].size() == 0){
                    level_block_components.push_back(block_components);
                    continue;
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
                    // std::cout << "prev_segments = " << +prev_segments << ", curr_segments = " << +curr_segments << std::endl; 
                    for(int i=prev_segments; i<curr_segments; i++){
                        // get id for blocks belonging to this aggregation block
                        const auto& global_blocks = aggregation_block_map.at(agg_block_id);
                        const uint8_t * compressed_data_pos = aggregated_level_segments[i - prev_segments];
                        // deal with each bitplane (which spans all blocks in the aggregation block)
                        // format: b1p1 b2p1 ... bnp1, b1p2 b2p2 ... bnp2, ...
                        // std::cout << "processing bitplanes in the segment\n";
                        // std::cout << "level = " << l << ", agg_block_id = " << b << std::endl;
                        for(int bp=0; bp<block_merge_counts[agg_block_id][i]; bp++){
                            uint8_t * precision_segment = NULL;
                            // printf("level = %d, agg_block = %d, segment = %d, bitplane = %d\n", l, agg_block_id, i, bp + bp_offset);
                            // printf("agg_block_bp_size = %d\n", agg_block_bp_sizes[agg_block_id][bp + bp_offset]);
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
                            // exit(0);
                            decompressed_segments.push_back(precision_segment);
                            // prepare each precision fragment
                            const uint8_t * bitplane_pos = precision_segment;
                            for(int j=0; j<global_blocks.size(); j++){
                                auto block_id = global_blocks[j];
                                auto block_num_element = block_level_elements[block_id][l];  
                                // TODO: change to adaptive computation with respect to T_stream
                                auto block_stream_size = compute_bitplane_size(block_num_element);                 
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
                    retriever.release();
                    // exit(0);                    
                }
                level_block_components.push_back(block_components);
            }
            total_retrieved_size += retrieved_size;
            printf("Retrieved size = %lu\n", retrieved_size);
            printf("Total retrieved size = %lu\n", total_retrieved_size);
            printf("Ratio = %.4f\n", num_elements * sizeof(T) * 1.0 / total_retrieved_size);
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
            deserialize(metadata_pos, num_dims, block_sizes);
            this->num_bitplanes = *(metadata_pos ++);
            num_levels = *(metadata_pos ++);
            std::vector<std::vector<int>> level_merge_counts;
            std::vector<std::vector<std::vector<double>>> level_squared_errors;
            std::vector<std::vector<uint32_t>> level_sizes;
            deserialize(metadata_pos, num_levels, level_max_exp);
            deserialize(metadata_pos, num_levels, level_merge_counts);
            deserialize(metadata_pos, num_levels, level_squared_errors);
            deserialize(metadata_pos, num_levels, level_aggregation_granularity);
            // deserialize(metadata_pos, num_levels, level_sizes);
            deserialize(metadata_pos, num_levels, level_num);
            deserialize(metadata_pos, num_levels, level_agg_block_bp_sizes);
            // read done
            num_blocks = std::vector<uint32_t>(dims.size());
            for(int i=0; i<dims.size(); i++){
                num_blocks[i] = (dims[i] - 1) / block_sizes[i] + 1;
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
            // compute level_block_merge_counts
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
            std::vector<uint32_t> block_dims(block_sizes);
            for(int i=0; i<nx; i++){
                block_dims[0] = (i < nx - 1) ? block_sizes[0] : (dims[0] - i*block_sizes[0]);
                for(int j=0; j<ny; j++){
                    block_dims[1] = (j < ny - 1) ? block_sizes[1] : (dims[1] - j*block_sizes[1]);
                    for(int k=0; k<nz; k++){
                        block_dims[2] = (k < nz - 1) ? block_sizes[2] : (dims[2] - k*block_sizes[2]);
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
            std::cout << "target_level = " << +target_level << std::endl;
            std::cout << "block_sizes = " << block_sizes[0] << " " << block_sizes[1] << " " << block_sizes[2] << " " << std::endl;
            // print_vec(num_blocks);
            auto nx = num_blocks[0];
            auto ny = num_blocks[1];
            auto nz = num_blocks[2];
            uint32_t total_num_blocks = nx * ny * nz;
            std::vector<uint32_t> block_dims(block_sizes);
            uint32_t max_num_block_elements = block_sizes[0] * block_sizes[1] * block_sizes[2];
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
                block_dims[0] = (i < nx - 1) ? block_sizes[0] : (dims[0] - i*block_sizes[0]);
                T * y_data_pos = x_data_pos;
                for(int j=0; j<ny; j++){
                    block_dims[1] = (j < ny - 1) ? block_sizes[1] : (dims[1] - j*block_sizes[1]);
                    T * z_data_pos = y_data_pos;
                    for(int k=0; k<nz; k++){
                        block_dims[2] = (k < nz - 1) ? block_sizes[2] : (dims[2] - k*block_sizes[2]);
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
                                    T * level_decoded_data = encoder.progressive_decode(components, level_elements[l], level_max_exp[l], prev_num_bp, curr_num_bp - prev_num_bp, l);
                                    const std::vector<uint32_t>& prev_dims = (l == 0) ? dims_dummy : level_dims[l - 1];
                                    interleaver.reposition(level_decoded_data, reconstruct_dimensions, level_dims[l], prev_dims, z_data_pos, this->strides);
                                    free(level_decoded_data);
                                }
                            }
                            // printf("recompose block %d\n", block_id);
                            // print_vec(block_dims);
                            // std::cout << reconstruct_dimensions[0] << " " << reconstruct_dimensions[1] << " " << reconstruct_dimensions[2] << std::endl;
                            // printf("target_level = %d\n", target_level);
                            decomposer.recompose(z_data_pos, reconstruct_dimensions, target_level, this->strides);
                        }
                        block_id ++;
                        z_data_pos += block_sizes[2];
                    }
                    y_data_pos += block_sizes[1] * dims[2];
                }
                x_data_pos += block_sizes[0] * dims[1] * dims[2];
            }
            return true;
        }

        inline uint32_t compute_bitplane_size(uint32_t block_num_element){
            return ((block_num_element - 1) / (UINT8_BITS*4) + 1) * sizeof(uint32_t);             
        }

        double compute_efficiency(int l, int i, int curr_num_segments, const std::vector<int>& blocks, double& error_gain, uint32_t& fetch_size){
            fetch_size = level_agg_block_sizes[l][i][curr_num_segments];
            // iterate through all blocks in the aggregation
            for(const auto& block_id : level_aggregation_block_map[l][i]){
                if(!converged[block_id]){
                    error_gain += error_estimator.estimate_error_gain(0, level_block_squared_errors[l][block_id][curr_num_segments], level_block_squared_errors[l][block_id][curr_num_segments+1], l);
                }
            }
            return error_gain / fetch_size;
        }

        Decomposer decomposer;
        Interleaver interleaver;
        Encoder encoder;
        SizeInterpreter interpreter;
        SNormErrorEstimator<T> error_estimator;
        Retriever retriever;
        Compressor compressor;
        std::vector<T> data;
        std::vector<uint32_t> dims;
        std::vector<int> level_aggregation_granularity;
        std::vector<int> level_max_exp;
        std::vector<std::vector<std::vector<int>>> level_block_merge_counts;            // L x sum(aggB)
        std::vector<std::vector<std::vector<double>>> level_block_squared_errors;       // L x n_b x P
        std::vector<std::vector<std::vector<const uint8_t*>>> level_block_components;   // L x n_b x P
        std::vector<std::vector<std::vector<uint32_t>>> level_agg_block_sizes;          // L x sum(aggB)
        std::vector<std::vector<std::vector<uint32_t>>> level_agg_block_bp_sizes;       // L x n_b x P
        std::vector<bool> converged;
        std::vector<uint32_t> level_num;
        std::vector<uint32_t> num_blocks;
        std::vector<std::vector<uint32_t>> level_num_agg_blocks;            // L x n_d
        std::vector<std::vector<uint8_t>> level_block_num_bitplanes;        // L x n_b
        std::vector<std::vector<uint8_t>> level_block_num_segments;         // sum(aggB)
        std::vector<std::vector<std::vector<uint32_t>>> block_level_dims;   // n_b x L x n_d
        std::vector<std::vector<uint32_t>> block_level_elements;            // n_b x L
        std::vector<uint32_t> strides;
        std::vector<uint32_t> block_sizes;
        uint8_t num_levels = 0;
        int num_bitplanes = 0;
        int num_elements = 0;
        int current_level = 0;
        size_t total_retrieved_size = 0;
        std::vector<std::unordered_map<int, AggregationBlockInfo>> level_block_aggregation_map; // map global block to aggregation block
        std::vector<std::unordered_map<int, std::vector<int>>> level_aggregation_block_map;     // map aggregation block to global blocks
    };
}
#endif

