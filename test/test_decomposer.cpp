#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include "utils.hpp"
#include "Decomposer/Decomposer.hpp"

using namespace std;

template <class T, class Decomposer>
void evaluate(const vector<T>& data, const vector<uint32_t>& dims, int target_level, Decomposer decomposer){
    struct timespec start, end;
    int err = 0;

    cout << "Using ";
    decomposer.print();
    cout << "Decompose data to " << target_level << " + 1 levels" << endl;
    vector<T> data_dup(data);
    err = clock_gettime(CLOCK_REALTIME, &start);
    decomposer.decompose(data_dup.data(), dims, target_level);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Decompose time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl; 

    err = clock_gettime(CLOCK_REALTIME, &start);
    decomposer.recompose(data_dup.data(), dims, target_level);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Recompose time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;

    T max_err = 0;
    for(int i=0; i<data.size(); i++){
        if(fabs(data[i] - data_dup[i]) > max_err){
            max_err = fabs(data[i] - data_dup[i]);
        }
    }
    cout << "Max error = " << max_err << endl;
}

template <class T>
void test(string filename, const vector<uint32_t>& dims){
    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    for(int target_level=0; target_level<5; target_level += 2){
        evaluate<T>(data, dims, target_level, MDR::MGARDOrthoganalDecomposer<T>());
        evaluate<T>(data, dims, target_level, MDR::MGARDHierarchicalDecomposer<T>());
    }
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