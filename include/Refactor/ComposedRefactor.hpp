#ifndef _MDR_COMPOSED_REFACTOR_HPP
#define _MDR_COMPOSED_REFACTOR_HPP

#include "RefactorInterface.hpp"
#include "Decomposer/Decomposer.hpp"
#include "Interleaver/Interleaver.hpp"
#include "BitplaneEncoder/BitplaneEncoder.hpp"
#include "ErrorCollector/ErrorCollector.hpp"
#include "ErrorEstimator/ErrorEstimator.hpp"
#include "LosslessCompressor/LevelCompressor.hpp"
#include "SizeInterpreter/SizeInterpreter.hpp"
#include "Writer/Writer.hpp"
#include "RefactorUtils.hpp"

namespace MDR {
    // a decomposition-based scientific data refactor: compose a refactor using decomposer, interleaver, encoder, and error collector
    template<class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorCollector, class Writer>
    class ComposedRefactor : public concepts::RefactorInterface<T> {
    public:
        ComposedRefactor(Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorCollector collector, Writer writer)
            : decomposer(decomposer), interleaver(interleaver), encoder(encoder), compressor(compressor), collector(collector), writer(writer) {}

        void refactor(T const * data_, const std::vector<uint32_t>& dims, uint8_t target_level, uint8_t num_bitplanes){
            dimensions = dims;
            uint32_t num_elements = 1;
            for(const auto& dim:dimensions){
                num_elements *= dim;
            }
            data = std::vector<T>(data_, data_ + num_elements);
            // if refactor successfully
            if(refactor(target_level, num_bitplanes)){
                // writer.write_level_components(level_components, level_sizes);
            }
            {
                level_abs_errors.clear();
                MaxErrorCollector<T> collector = MaxErrorCollector<T>();
                for(int i=0; i<=target_level; i++){
                    auto collected_error = collector.collect_level_error(NULL, 0, level_squared_errors[i].size(), level_error_bounds[i]);
                    level_abs_errors.push_back(collected_error);
                }                
            }
        }

        uint8_t * write_metadata(uint32_t& size) const {
            uint32_t metadata_size = sizeof(uint8_t) + get_size(dimensions) // dimensions
                            + sizeof(uint8_t) + get_size(level_error_bounds) + get_size(level_squared_errors) + get_size(level_sizes); // level information
            uint8_t * metadata = (uint8_t *) malloc(metadata_size);
            uint8_t * metadata_pos = metadata;
            *(metadata_pos ++) = (uint8_t) dimensions.size();
            serialize(dimensions, metadata_pos);
            *(metadata_pos ++) = (uint8_t) level_error_bounds.size();
            serialize(level_error_bounds, metadata_pos);
            serialize(level_squared_errors, metadata_pos);
            serialize(level_sizes, metadata_pos);
            // writer.write_metadata(metadata, metadata_size);
            size = metadata_size;
            return metadata;
        }

        uint8_t * get_data(T value_range, std::vector<int>& positions, uint32_t& size) {
            uint8_t * reordered_data = reorganize(value_range, positions, size);
            return reordered_data;
        }

        ~ComposedRefactor(){
            for(int i=0; i<level_components.size(); i++){
                for(int j=0; j<level_components[i].size(); j++){
                    free(level_components[i][j]);                    
                }
            }            
        }

        void print() const {
            std::cout << "Composed refactor with the following components." << std::endl;
            std::cout << "Decomposer: "; decomposer.print();
            std::cout << "Interleaver: "; interleaver.print();
            std::cout << "Encoder: "; encoder.print();
        }
    private:
        uint8_t * reorganize(T value_range, std::vector<int>& positions, uint32_t& total_size){
            positions.clear();
            std::vector<double> eb{0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0};
            for(int i=0; i<eb.size(); i++){
                eb[i] *= value_range;
            }
            int eb_index = 0;
            auto error_estimator = MDR::MaxErrorEstimatorHB<T>();
            const int num_levels = level_sizes.size();
            total_size = 0;
            for(int i=0; i<num_levels; i++){
                for(int j=0; j<level_sizes[i].size(); j++){
                    total_size += level_sizes[i][j];
                }
            }
            const std::vector<std::vector<double>>& level_errors = level_abs_errors;
            uint8_t * reorganized_data = (uint8_t *) malloc(total_size);
            uint8_t * reorganized_data_pos = reorganized_data;
            std::vector<int> index(level_sizes.size(), 0);
            double accumulated_error = 0;
            std::priority_queue<UnitErrorGain, std::vector<UnitErrorGain>, CompareUniteErrorGain> heap;
            for(int i=0; i<num_levels; i++){
                memcpy(reorganized_data_pos, level_components[i][0], level_sizes[i][0]);
                reorganized_data_pos += level_sizes[i][0];
                index[i] ++;
                accumulated_error += error_estimator.estimate_error(level_errors[i][index[i]], i);
                if(accumulated_error < eb[eb_index]){
                    positions.push_back(reorganized_data_pos - reorganized_data);
                    eb_index ++;
                }
                std::cout << i;
                // push the next one
                if(index[i] != level_sizes[i].size()){
                    double error_gain = error_estimator.estimate_error_gain(0, level_errors[i][index[i]], level_errors[i][index[i] + 1], i);
                    heap.push(UnitErrorGain(error_gain / level_sizes[i][index[i]], i));
                }
            }
            while(!heap.empty()){
                auto unit_error_gain = heap.top();
                heap.pop();
                int i = unit_error_gain.level;
                int j = index[i];
                memcpy(reorganized_data_pos, level_components[i][j], level_sizes[i][j]);
                reorganized_data_pos += level_sizes[i][j];
                index[i] ++;
                accumulated_error -= error_estimator.estimate_error(level_errors[i][j], i);
                accumulated_error += error_estimator.estimate_error(level_errors[i][j+1], i);
                if(accumulated_error < eb[eb_index]){
                    positions.push_back(reorganized_data_pos - reorganized_data);
                    eb_index ++;
                }
                if(index[i] != level_sizes[i].size()){
                    double error_gain = error_estimator.estimate_error_gain(0, level_errors[i][index[i]], level_errors[i][index[i] + 1], i);
                    heap.push(UnitErrorGain(error_gain / level_sizes[i][index[i]], i));
                }
                std::cout << i;
            }
            std::cout << std::endl;
            return reorganized_data;
        }
        bool refactor(uint8_t target_level, uint8_t num_bitplanes){
            uint8_t max_level = log2(*min_element(dimensions.begin(), dimensions.end())) - 1;
            if(target_level > max_level){
                std::cerr << "Target level is higher than " << max_level << std::endl;
                return false;
            }
            // decompose data hierarchically
            decomposer.decompose(data.data(), dimensions, target_level);

            // encode level by level
            level_error_bounds.clear();
            level_squared_errors.clear();
            level_components.clear();
            level_sizes.clear();
            auto level_dims = compute_level_dims(dimensions, target_level);
            auto level_elements = compute_level_elements(level_dims, target_level);
            std::vector<uint32_t> dims_dummy(dimensions.size(), 0);
            SquaredErrorCollector<T> s_collector = SquaredErrorCollector<T>();
            for(int i=0; i<=target_level; i++){
                const std::vector<uint32_t>& prev_dims = (i == 0) ? dims_dummy : level_dims[i - 1];
                T * buffer = (T *) malloc(level_elements[i] * sizeof(T));
                // extract level i component
                interleaver.interleave(data.data(), dimensions, level_dims[i], prev_dims, reinterpret_cast<T*>(buffer));
                // compute max coefficient as level error bound
                T level_max_error = compute_max_abs_value(reinterpret_cast<T*>(buffer), level_elements[i]);
                level_error_bounds.push_back(level_max_error);
                // collect errors
                auto collected_error = s_collector.collect_level_error(buffer, level_elements[i], num_bitplanes, level_max_error);
                level_squared_errors.push_back(collected_error);
                // encode level data
                int level_exp = 0;
                frexp(level_max_error, &level_exp);
                std::vector<uint32_t> stream_sizes;
                auto streams = encoder.encode(buffer, level_elements[i], level_exp, num_bitplanes, stream_sizes);
                free(buffer);
                // lossless compression
                compressor.compress_level(streams, stream_sizes);
                // record encoded level data and size
                level_components.push_back(streams);
                level_sizes.push_back(stream_sizes);
            }
            return true;
        }

        Decomposer decomposer;
        Interleaver interleaver;
        Encoder encoder;
        Compressor compressor;
        ErrorCollector collector;
        Writer writer;
        std::vector<T> data;
        std::vector<uint32_t> dimensions;
        std::vector<T> level_error_bounds;
        std::vector<int> order;
        std::vector<std::vector<uint8_t*>> level_components;
        std::vector<std::vector<uint32_t>> level_sizes;
        std::vector<std::vector<double>> level_squared_errors;
        std::vector<std::vector<double>> level_abs_errors;
    };
}
#endif

