#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include "utils.hpp"
#include "BitplaneEncoder/BitplaneEncoder.hpp"
#include "ErrorCollector/MaxErrorCollector.hpp"
#include "ErrorCollector/SNormErrorCollector.hpp"
#include "ErrorEstimator/MaxErrorEstimator.hpp"
#include "ErrorEstimator/SNormErrorEstimator.hpp"
#include "Retriever/GeneralRetriever.hpp"

#define MAX(a, b) (a>b) ? (a) : (b)
using namespace std;

template <class T, class ErrorEstimator>
void evaluate(const vector<T>& data, int level_exp, const vector<const uint8_t*>& streams, vector<vector<uint32_t>>& level_sizes, vector<vector<double>>& level_errors, ErrorEstimator estimator){
    struct timespec start, end;
    int err = 0;

    auto retriever = MDR::GeneralRetriever<ErrorEstimator>(estimator);
    cout << "Using ";
    estimator.print();
    cout << "Using ";
    retriever.print();

    cout << endl;
    cout << "level_errors: " << endl;
    for(int i=0; i<level_errors.size(); i++){
        for(int j=0; j<level_errors[i].size(); j++){
            cout << level_errors[i][j] << " ";
        }
        cout << endl;
    }
    uint32_t num_elements = data.size();
    vector<uint8_t> order(level_sizes[0].size());
    const int test_num = 4;
    double target_errors[test_num] = {1.0, 1e-2, 1e-4, 1e-6};
    MDR::GroupedBPEncoder<T, uint32_t> encoder;
    for(int i=0; i<test_num; i++){
        cout << "Test " << i << ", target_error = " <<target_errors[i] << endl;
        fflush(stdout);
        vector<uint8_t> index;
        err = clock_gettime(CLOCK_REALTIME, &start);
        uint32_t retrieve_size = retriever.interpret_size(level_sizes, level_errors, order, target_errors[i], index);
        err = clock_gettime(CLOCK_REALTIME, &end);        
        cout << "Retrieve size = " << retrieve_size << ", level 0 #bitplane = " << +index[0] << endl;
        cout << "Interpret size time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
        auto dec_data = encoder.decode(streams, num_elements, level_exp, index[0]);
        double max_error = 0;
        double squared_error = 0;
        for(int i=0; i<num_elements; i++){
            max_error = MAX(max_error, fabs(data[i] - dec_data[i]));
            squared_error += (data[i] - dec_data[i]) * (data[i] - dec_data[i]);
        }
        free(dec_data);
        cout << "Max error = " << max_error << ", squared error = " << squared_error << endl;
    }
    cout << endl;
}

template <class T>
void test(string filename){
    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    T max_value = 0;
    for(int i=0; i<num_elements; i++){
        if(fabs(data[i]) > max_value) max_value = fabs(data[i]);
    }
    int level_exp = 0;
    frexp(max_value, &level_exp);

    const int num_bitplanes = 32;
    vector<vector<uint32_t>> level_sizes;
    vector<uint32_t> sizes;
    MDR::GroupedBPEncoder<T, uint32_t> encoder;
    vector<uint8_t*> streams = encoder.encode(data.data(), num_elements, level_exp, num_bitplanes, sizes);
    level_sizes.push_back(sizes);

    std::vector<uint8_t const*> streams_const;
    for(int i=0; i<streams.size(); i++){
        streams_const.push_back(streams[i]);
    }

    vector<vector<double>> level_errors;
    MDR::MaxErrorCollector<T> collector;
    auto errors = collector.collect_level_error(data.data(), num_elements, num_bitplanes, max_value);
    level_errors.push_back(errors);
    vector<vector<double>> level_s_errors;
    MDR::SNormErrorCollector<T> collector_s;
    auto s_errors = collector_s.collect_level_error(data.data(), num_elements, num_bitplanes, max_value);
    level_s_errors.push_back(s_errors);

    MDR::MaxErrorEstimatorHB<T> max_e_estimator = MDR::MaxErrorEstimatorHB<T>();
    MDR::SNormErrorEstimator<T> s_norm_e_estimator = MDR::SNormErrorEstimator<T>(1, 0, 0);
    // evaluate(data, level_exp, streams_const, level_sizes, level_errors, max_e_estimator);
    evaluate(data, level_exp, streams_const, level_sizes, level_s_errors, s_norm_e_estimator);
    for(int i=0; i<streams.size(); i++){
        free(streams[i]);
    }
}

int main(int argc, char ** argv){

    string filename = string(argv[1]);
    test<float>(filename);
    return 0;

}