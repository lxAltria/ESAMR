#ifndef _MDR_ADAPTIVE_REFACTOR_HPP
#define _MDR_ADAPTIVE_REFACTOR_HPP

#include "RefactorInterface.hpp"
#include "Decomposer/Decomposer.hpp"
#include "Interleaver/Interleaver.hpp"
#include "BitplaneEncoder/BitplaneEncoder.hpp"
#include "ErrorCollector/ErrorCollector.hpp"
#include "LosslessCompressor/ZSTD.hpp"
#include "Writer/Writer.hpp"
#include "RefactorUtils.hpp"

namespace MDR {
    // a decomposition-based scientific data refactor: compose a refactor using decomposer, interleaver, encoder, and error collector
    template<class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorCollector, class Writer>
    class AdaptiveRefactor : public concepts::RefactorInterface<T> {
    public:
        AdaptiveRefactor(Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorCollector collector, Writer writer)
            : decomposer(decomposer), interleaver(interleaver), encoder(encoder), compressor(compressor), collector(collector), writer(writer) {}

        void refactor(T const * data_, const std::vector<uint32_t>& dims, uint8_t target_level, uint8_t num_bitplanes, uint32_t block_size=0){
            if(block_size == 0){
                block_size = dims[0];
                for(int i=1; i<dims.size(); i++){
                    if(block_size > dims[i]) block_size = dims[i];
                }
            }
            this->block_size = block_size;
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

            this->dims = dims;
            uint32_t num_elements = 1;
            for(const auto& dim:dims){
                num_elements *= dim;
            }
            data = std::vector<T>(data_, data_ + num_elements);
            // if refactor successfully
            refactor(target_level, num_bitplanes);
        }

        void write_metadata() const {
        }

        ~AdaptiveRefactor(){}

        void print() const {
            std::cout << "Adaptive refactor with the following components." << std::endl;
            std::cout << "Decomposer: "; decomposer.print();
            std::cout << "Interleaver: "; interleaver.print();
            std::cout << "Encoder: "; encoder.print();
        }
    private:
        // global dims and strides known as class variable
        void decompose_and_interleave(T * data, int target_level, const std::vector<uint32_t>& block_dims, const std::vector<std::vector<uint32_t>>& level_dims, 
                                        const std::vector<uint32_t>& level_elements, T * buffer_pos){
            const std::vector<uint32_t> dims_dummy(dims.size(), 0);
            decomposer.decompose(data, block_dims, target_level, this->strides);
            // std::cout << "decompose done" << std::endl;
            // print_vec(level_elements);
            T * cur_buffer_pos = buffer_pos;
            for(int l=0; l<=target_level; l++){
                const std::vector<uint32_t>& prev_dims = (l == 0) ? dims_dummy : level_dims[l - 1];
                // use dimensions instead of block_dims because of global strides
                interleaver.interleave(data, this->dims, level_dims[l], prev_dims, reinterpret_cast<T*>(cur_buffer_pos), this->strides);
                cur_buffer_pos += level_elements[l];
            }            
            // std::cout << "interleave done" << std::endl;
        }

        // T compute_max_level_value(T const * level_data, int num_blocks, ptrdiff_t block_offset, const std::vector<std::vector<uint32_t>>& block_level_elements){
        //     T max_val = 0;
        //     T const * level_data_pos = level_data;
        //     for(int i=0; i<num_blocks; i++){
        //         for(int j=0; j<block_level_elements[i]; j++){
        //             if(max_val < fabs(level_data_pos[j])) max_val = level_data_pos[j];
        //         }
        //         level_data_pos += block_offset;
        //     }
        //     return max_val;
        // }

        void collect_components(int target_level, int aggregation_size, uint32_t max_num_block_elements, uint32_t num_bitplanes,
                                    const std::vector<std::vector<uint32_t>>& block_level_elements, const std::vector<uint32_t>& level_offset, T const * data){
            std::cout << "collect_components start" << std::endl;
            int aggregation_total_num = aggregation_size * aggregation_size * aggregation_size;
            // buffer is used to store bitplanes
            uint8_t * buffer = (uint8_t *) malloc(SEGMENT_SIZE * 1024 * sizeof(T));
            std::vector<int> bitplane_sizes(num_blocks[0]*num_blocks[1]*num_blocks[2], 0);
            int segment_count = 0;
            print_vec(level_offset);
            // aggregation_granularity indicates how many subregions are used along each dimension
            // start at finest granularity
            // i.e., all blocks are in the same aggregation area
            int aggregation_granularity = 1;
            for(int l=0; l<=target_level; l++){
                std::cout << "level " << l << std::endl;
                std::vector<uint8_t*> level_component;
                std::vector<uint32_t> level_size;
                // compute number of aggregations based on aggregation granularity
                // nx, ny, nz: number of aggregation blocks along each dimension
                int nx = (num_blocks[0] < aggregation_granularity) ? num_blocks[0] : aggregation_granularity;
                int ny = (num_blocks[1] < aggregation_granularity) ? num_blocks[1] : aggregation_granularity;
                int nz = (num_blocks[2] < aggregation_granularity) ? num_blocks[2] : aggregation_granularity;
                // agg_block_size_x, agg_block_size_y, agg_block_size_z: number of data blocks in each aggregation block
                int agg_block_size_x = num_blocks[0] / nx;
                int agg_block_size_y = num_blocks[1] / ny;
                int agg_block_size_z = num_blocks[2] / nz;
                // residual for computing block size
                int residual_x = num_blocks[0] % nx;
                int residual_y = num_blocks[1] % ny;
                int residual_z = num_blocks[2] % nz;
                std::cout << "aggregation dimensions: ";
                std::cout << nx << " " << ny << " " << nz << std::endl;
                std::cout << "aggregation block_size: ";
                std::cout << agg_block_size_x << " " << agg_block_size_y << " " << agg_block_size_z << std::endl;
                // iterate each aggregation area
                for(int i=0; i<nx; i++){
                    int aggregation_nx = (i < residual_x) ? (agg_block_size_x + 1) : agg_block_size_x;
                    int block_offset_x = (i < residual_x) ? (i * (agg_block_size_x + 1)) : (i*agg_block_size_x + residual_x);
                    for(int j=0; j<ny; j++){
                        int aggregation_ny = (j < residual_y) ? (agg_block_size_y + 1) : agg_block_size_y;
                        int block_offset_y = (j < residual_y) ? (j * (agg_block_size_y + 1)) : (i*agg_block_size_y + residual_y);
                        for(int k=0; k<nz; k++){
                            int aggregation_nz = (k < residual_z) ? (agg_block_size_z + 1) : agg_block_size_z;
                            int block_offset_z = (k < residual_z) ? (k * (agg_block_size_z + 1)) : (k*agg_block_size_z + residual_z);
                            int num_aggregation = aggregation_nx * aggregation_ny * aggregation_nz;
                            // iterate data in one aggregation area which has aggregation_size^3 blocks
                            // buffer is used to store the encoded data
                            std::cout << "aggregation num_blocks: ";
                            std::cout << aggregation_nx << " " << aggregation_ny << " " << aggregation_nz << std::endl;
                            uint8_t * bitplane_pos = buffer;
                            int aggregation_id = 0;
                            for(int ii=0; ii<aggregation_nx; ii++){
                                for(int jj=0; jj<aggregation_ny; jj++){
                                    for(int kk=0; kk<aggregation_nz; kk++){
                                        // deal with data in one block 
                                        // identify the location of data: block_offset + level_offset
                                        int block_id = (block_offset_x + ii) * num_blocks[1] * num_blocks[2] + (block_offset_y + jj) * num_blocks[2] + (block_offset_z + kk);
                                        std::cout << "block_id = " << block_id << "\n";
                                        auto num_level_elements = block_level_elements[block_id][l];
                                        T const * data_pos = data + block_id * max_num_block_elements + level_offset[l];
                                        // use negabinary encoder, where each bitplane has the same stream_size
                                        // TODO: consider other encoders
                                        // TODO: use a global level max error
                                        T level_max_error = compute_max_abs_value(reinterpret_cast<T const*>(data_pos), num_level_elements);
                                        int level_exp = 0;
                                        frexp(level_max_error, &level_exp);
                                        std::vector<uint32_t> stream_sizes;
                                        std::vector<double> level_sq_err;
                                        std::cout << "start encoding: " << "offset = " << data_pos - data << ", num_level_elements = " << num_level_elements << "\n";
                                        auto streams = encoder.encode(data_pos, num_level_elements, level_exp, num_bitplanes, stream_sizes, level_sq_err);
                                        bitplane_sizes[aggregation_id] = stream_sizes[0];
                                        std::cout << "bitplane_size = " << stream_sizes[0] << "\n";
                                        // copy data for alignment on bitplanes
                                        uint8_t * current_bitplane_pos = bitplane_pos;
                                        for(int bp=0; bp<num_bitplanes; bp++){
                                            memcpy(current_bitplane_pos, streams[bp], stream_sizes[bp]);
                                            current_bitplane_pos += bitplane_sizes[aggregation_id] * num_aggregation;
                                            free(streams[bp]);
                                        }
                                        bitplane_pos += bitplane_sizes[aggregation_id];
                                        aggregation_id ++;
                                    }
                                }
                            }
                            std::cout << std::endl;
                            // TODO: bitplane size is larger than SEGMENT_SIZE 
                            uint8_t * level_component_buffer = (uint8_t *) malloc(SEGMENT_SIZE);
                            // compress each bitplane
                            int aggregated_bitplane_size = 0;
                            for(int i=0; i<num_aggregation; i++){
                                aggregated_bitplane_size += bitplane_sizes[i];
                            }
                            uint32_t current_size = 0;
                            uint8_t * buffer_pos = buffer;
                            uint8_t * level_component_buffer_pos = level_component_buffer;
                            std::cout << "aggregated_bitplane_size = " << aggregated_bitplane_size << std::endl;
                            for(int bp=0; bp<num_bitplanes; bp++){
                                uint8_t * compressed_data = NULL;
                                auto compressed_size = ZSTD::compress(buffer_pos, aggregated_bitplane_size, &compressed_data);
                                // concatenate bitplane and put as components
                                if(current_size + compressed_size > SEGMENT_SIZE){
                                    level_component.push_back(level_component_buffer);
                                    level_size.push_back(current_size);
                                    level_component_buffer = (uint8_t *) malloc(SEGMENT_SIZE);
                                    current_size = 0;
                                    level_component_buffer_pos = level_component_buffer;
                                    segment_count ++;                                  
                                }
                                memcpy(level_component_buffer_pos, compressed_data, compressed_size);
                                free(compressed_data);
                                current_size += compressed_size;
                                buffer_pos += aggregated_bitplane_size;
                            }
                            level_component.push_back(level_component_buffer);
                            level_size.push_back(current_size);
                        }
                    }
                }
                level_components.push_back(level_component);
                level_sizes.push_back(level_size);
                // check whether to change aggregation_size
                if(segment_count){
                    aggregation_granularity *= 2;
                }
            }
            // TODO: free level_components in other place
            for(int i=0; i<level_components.size(); i++){
                free(level_components[i]);
            }
            free(buffer);
        }

        bool refactor(uint8_t target_level, uint8_t num_bitplanes){
            uint8_t max_level = log2(*min_element(dims.begin(), dims.end())) - 1;
            if(target_level > max_level){
                std::cerr << "Target level is higher than " << max_level << std::endl;
                return false;
            }
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
            std::vector<std::vector<uint32_t>> block_level_elements;
            // decompose and interleave
            // the buffer store data in a blockwise then levelwise fashtion
            // i.e.: |B1L1, B1L2, ... B2L1, B2L2 ...|
            // stride between different block is max_num_block_elements
            T * buffer = (T *) malloc(max_num_block_elements * total_num_blocks * sizeof(T));
            std::cout << "max_num_block_elements = " << max_num_block_elements << ", total_num_blocks = " << total_num_blocks << std::endl;
            T * buffer_pos = buffer;
            // int block_id = 0;
            T * x_data_pos = data.data();
            for(int i=0; i<nx; i++){
                block_dims[0] = (i < nx - 1) ? block_size : (dims[0] - i*block_size);
                T * y_data_pos = x_data_pos;
                for(int j=0; j<ny; j++){
                    block_dims[1] = (j < ny - 1) ? block_size : (dims[1] - j*block_size);
                    T * z_data_pos = y_data_pos;
                    for(int k=0; k<nz; k++){
                        block_dims[2] = (k < nz - 1) ? block_size : (dims[2] - k*block_size);
                        // print_vec(block_dims);
                        // decompose and interleave data in each block
                        decompose_and_interleave(z_data_pos, target_level, block_dims, max_level_dims, max_level_elements, buffer_pos);
                        // std::cout << "block decompose and interleave done" << std::endl;
                        auto level_dims = compute_level_dims(block_dims, target_level);
                        auto level_elements = compute_level_elements(level_dims, target_level);
                        block_level_elements.push_back(level_elements);
                        // block_id ++;
                        buffer_pos += max_num_block_elements;
                        z_data_pos += block_size;
                    }
                    y_data_pos += block_size * dims[2];
                }
                x_data_pos += block_size * dims[1] * dims[2];
            }
            std::cout << "decompose and interleave done" << std::endl;
            std::vector<uint32_t> level_offset;
            uint32_t offset = 0;
            for(int l=0; l<=target_level; l++){
                level_offset.push_back(offset);
                offset += max_level_elements[l];
            }
            // very coarse granularity for the first three levels
            // use aggregated exponent for each level
            // use adaptive granularity
            int aggregation_size = 1;
            collect_components(target_level, aggregation_size, max_num_block_elements, num_bitplanes, block_level_elements, level_offset, buffer);
            std::cout << "collect_components done" << std::endl;
            free(buffer);
            return true;
        }

        Decomposer decomposer;
        Interleaver interleaver;
        Encoder encoder;
        Compressor compressor;
        ErrorCollector collector;
        Writer writer;
        std::vector<T> data;
        std::vector<uint32_t> dims;
        std::vector<T> level_error_bounds;
        std::vector<uint8_t> stopping_indices;
        std::vector<std::vector<uint8_t*>> level_components;
        std::vector<std::vector<uint32_t>> level_sizes;
        std::vector<uint32_t> level_num;
        std::vector<std::vector<double>> level_squared_errors;
        const uint32_t SEGMENT_SIZE = 1024 * 1024; // 1 MB
        std::vector<uint32_t> num_blocks;
        std::vector<uint32_t> strides;
        int block_size;
    };
}
#endif

