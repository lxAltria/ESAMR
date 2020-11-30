#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include "utils.hpp"
#include "BitplaneEncoder/BitplaneEncoder.hpp"
#include "ErrorCollector/ErrorCollector.hpp"

#define MAX(a, b) (a>b) ? (a) : (b)
using namespace std;

template <class T>
void compute_true_errors(const vector<T>& data, int num_elements, T max_value){
    const int num_bitplanes = 32;
    vector<double> true_max_error(num_bitplanes + 1, 0);
    vector<double> true_squared_error(num_bitplanes + 1, 0);
    MDR::GroupedBPEncoder<T, uint32_t> encoder;
    int level_exp = 0;
    frexp(max_value, &level_exp);
    vector<uint32_t> sizes;
    vector<uint8_t*> streams = encoder.encode(data.data(), num_elements, level_exp, num_bitplanes, sizes);
    std::vector<uint8_t const*> streams_const;
    for(int i=0; i<streams.size(); i++){
        streams_const.push_back(streams[i]);
    }
    for(int i=1; i<=num_bitplanes; i++){
        auto dec_data = encoder.decode(streams_const, num_elements, level_exp, i);
        double max_error = 0;
        double squared_error = 0;
        for(int j=0; j<data.size(); j++){
            max_error = MAX(max_error, fabs(data[j] - dec_data[j]));
            squared_error += (data[j] - dec_data[j]) * (data[j] - dec_data[j]);
        }
        true_max_error[i] = max_error;
        true_squared_error[i] = squared_error;
        // cout << "max_error = " << max_error << ", squared_error = " << squared_error << endl;
        free(dec_data);
    }
    for(int i=0; i<streams.size(); i++){
        free(streams[i]);
    }
    cout << "True max errors: " << endl;
    for(int i=0; i<true_max_error.size(); i++){
        cout << true_max_error[i] << " ";
    }
    cout << endl;

    cout << "True squared errors: " << endl;
    for(int i=0; i<true_squared_error.size(); i++){
        cout << true_squared_error[i] << " ";
    }
    cout << endl;
    cout << endl;
}

template <class T, class ErrorCollector>
void evaluate(const vector<T>& data, int num_elements, T max_value, ErrorCollector collector){
    struct timespec start, end;
    int err = 0;

    const int num_bitplanes = 32;
    cout << "Using ";
    collector.print();
    err = clock_gettime(CLOCK_REALTIME, &start);
    auto collected_error = collector.collect_level_error(data.data(), num_elements, num_bitplanes, max_value);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Error collection time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;

    cout << "Error vector size = " << collected_error.size() << endl;
    for(int i=0; i<collected_error.size(); i++){
        cout << collected_error[i] << " ";
    }
    cout << endl;
    cout << endl;
}

template <class T>
void test(string filename){
    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    T max_val = 0;
    for(int i=0; i<num_elements; i++){
        if(fabs(data[i]) > max_val) max_val = fabs(data[i]);
    }
    compute_true_errors(data, num_elements, max_val);
    evaluate(data, num_elements, max_val, MDR::MaxErrorCollector<T>());
    evaluate(data, num_elements, max_val, MDR::SquaredErrorCollector<T>());
}

int main(int argc, char ** argv){

    string filename = string(argv[1]);
    test<float>(filename);
    return 0;

}