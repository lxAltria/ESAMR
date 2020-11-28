#ifndef _MDR_GENERAL_RETRIEVER_HPP
#define _MDR_GENERAL_RETRIEVER_HPP

#include "RetrieverInterface.hpp"

namespace MDR {
    // General error-controlled data retriever
    template<class ErrorEstimator>
    class GeneralRetriever : public concepts::RetrieverInterface {
    public:
        GeneralRetriever(const ErrorEstimator& e){
            error_estimator = e;
        }
        uint32_t interpret_size(const std::vector<std::vector<uint32_t>>& level_sizes, const std::vector<std::vector<double>>& level_errors, const std::vector<uint8_t>& order, double tolerance, std::vector<uint8_t>& index) const {
            // #level_error[i] equals #bitplane + 1 as the first column records the error when reconstructed data are all 0
            size_t retrieve_size = 0;
            int count = 0;
            // note roundoff error for accumulated errors
            double accumulated_error = 0;
            uint8_t num_levels = level_errors.size();
            for(int i=0; i<num_levels; i++){
                accumulated_error += error_estimator.estimate_error(level_errors[i][0], i);
            }
            // for(int i=0; i<num_levels; i++){
            //     for(int j=0; j<level_errors[j].size(); j++){
            //         std::cout << level_errors[i][j] << " ";
            //     }
            //     std::cout << std::endl;
            // }
            // std::cout << std::endl;
            // for(int i=0; i<order.size(); i++){
            //     std::cout << +order[i] << " ";
            // }
            // std::cout << std::endl;
            index.clear();
            index.resize(num_levels, 0);
            while((accumulated_error > tolerance) && (count < order.size())){
                uint8_t level = order[count ++];
                uint8_t bitplane_index = index[level];
                retrieve_size += level_sizes[level][bitplane_index];
                accumulated_error -= error_estimator.estimate_error(level_errors[level][bitplane_index], level);
                accumulated_error += error_estimator.estimate_error(level_errors[level][bitplane_index + 1], level);
                index[level] ++;
            }
            std::cout << "Retrieve_size = " << retrieve_size << std::endl;
            std::cout << "Estimate error = " << accumulated_error << ", tolerance = " << tolerance << std::endl;
            return retrieve_size;
        }
        void print() const {
            std::cout << "General retriever." << std::endl;
        }
    private:
        ErrorEstimator error_estimator;
    };
}
#endif
