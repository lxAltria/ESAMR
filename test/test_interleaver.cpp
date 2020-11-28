#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include "utils.hpp"
#include "RefactorUtils.hpp"
#include "Interleaver/DirectInterleaver.hpp"

using namespace std;

template <class T, class Interleaver>
void evaluate(const vector<T>& data, const vector<uint32_t>& dims, int target_level, Interleaver interleaver){
    struct timespec start, end;
    int err = 0;
    cout << "Using ";
    interleaver.print();
    auto level_dims = MDR::compute_level_dims(dims, target_level);
    auto level_elements = MDR::compute_level_elements(level_dims, target_level);
    vector<uint32_t> dims_dummy(dims.size(), 0);
    vector<uint8_t *> buffers;
    for(int i=0; i<=target_level; i++){
        err = clock_gettime(CLOCK_REALTIME, &start);
        const vector<uint32_t>& prev_dims = (i == 0) ? dims_dummy : level_dims[i - 1];
        uint8_t * buffer = (uint8_t *) malloc(level_elements[i] * sizeof(T));
        interleaver.interleave(data.data(), dims, level_dims[i], prev_dims, reinterpret_cast<T*>(buffer));
        buffers.push_back(buffer);
        err = clock_gettime(CLOCK_REALTIME, &end);
        cout << "Level " << target_level - i << " interleave time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl; 
    }
    vector<uint8_t const *> buffers_const;
    for(int i=0; i<=target_level; i++){
        buffers_const.push_back(buffers[i]);
    }
    vector<T> data_reposition(data.size(), 0);
    for(int i=0; i<=target_level; i++){
        err = clock_gettime(CLOCK_REALTIME, &start);
        const vector<uint32_t>& prev_dims = (i == 0) ? dims_dummy : level_dims[i - 1];
        interleaver.reposition(reinterpret_cast<const T *>(buffers_const[i]), dims, level_dims[i], prev_dims, data_reposition.data());
        err = clock_gettime(CLOCK_REALTIME, &end);
        cout << "Level " << target_level - i << " reposition time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    }
    for(int i=0; i<buffers.size(); i++){
        free(buffers[i]);
    }
    T max_err = 0;
    for(int i=0; i<data.size(); i++){
        if(fabs(data[i] - data_reposition[i]) > max_err){
            max_err = fabs(data[i] - data_reposition[i]);
        }
    }
    cout << "Max error = " << max_err << endl;
}

template <class T>
void test(string filename, const vector<uint32_t>& dims){
    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    const int target_level = 3;
    evaluate<T>(data, dims, target_level, MDR::DirectInterleaver<T>());
}

int main(int argc, char ** argv){

    string filename = string(argv[1]);
    int num_dims = atoi(argv[2]);
    vector<uint32_t> dims(num_dims);
    for(int i=0; i<num_dims; i++){
        dims[i] = atoi(argv[3+i]);
    }
    test<float>(filename, dims);
    return 0;

}