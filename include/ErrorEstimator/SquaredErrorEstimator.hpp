#ifndef _MDR_SQUARED_ERROR_ESTIMATOR_HPP
#define _MDR_SQUARED_ERROR_ESTIMATOR_HPP

#include "ErrorEstimatorInterface.hpp"

namespace MDR {
    template<class T>
    class SquaredErrorEstimator : public concepts::ErrorEstimatorInterface<T>{};

    // S-norm error estimator for orthogonal basis
    template<class T>
    class SNormErrorEstimator : public SquaredErrorEstimator<T> {
    public:
        SNormErrorEstimator(int num_dims, int target_level, T s) : s(s) {
            s_table = std::vector<T>(target_level + 1);
            for(int i=0; i<=target_level; i++){
                // 2^(sl) * vol(P_l) where vol(P_l) = 2^(dl)
                int l = target_level - i;
                s_table[i] = pow(2, 2*s*l + num_dims*l);
            }
        }
        SNormErrorEstimator() : SNormErrorEstimator(1, 0, 0) {}
        inline T estimate_error(T error, int level) const {
            return s_table[level] * error;
        }
        inline T estimate_error(T data, T reconstructed_data, int level) const {
            return s_table[level] * (data - reconstructed_data);
        }
        inline T estimate_error_gain(T base, T current_level_err, T next_level_err, int level) const {
            return s_table[level] * (current_level_err - next_level_err);
        }
        void print() const {
            std::cout << "S-norm error estimator (s = " << s << ")." << std::endl;
        }
    private:
        T s = 0;
        std::vector<T> s_table;
    };
}
#endif
