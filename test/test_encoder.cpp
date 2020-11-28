#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include "utils.hpp"
#include "BitplaneEncoder/BitplaneEncoder.hpp"

using namespace std;

template <class T, class T2>
void evaluate(const vector<T>& data, size_t num_elements, int level_exp, int num_bitplanes){
    struct timespec start, end;
    int err = 0;

    MDR::GroupedBPEncoder<T, T2> encoder;
    vector<uint32_t> sizes;
    err = clock_gettime(CLOCK_REALTIME, &start);
    std::vector<uint8_t*> streams = encoder.encode(data.data(), num_elements, level_exp, num_bitplanes, sizes);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << streams.size() << endl;
    cout << "Encoding time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;	

    std::vector<uint8_t const*> streams_const;
    for(int i=0; i<streams.size(); i++){
        streams_const.push_back(streams[i]);
    }
    err = clock_gettime(CLOCK_REALTIME, &start);
    auto dec_data = encoder.decode(streams_const, num_elements, level_exp, num_bitplanes);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Decoding time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;

    cout << "Encoded sizes: ";
    for(int i=0; i<sizes.size(); i++){
    	cout << sizes[i] << " ";
        free(streams[i]);
    }
    cout << endl;

    T max_err = 0;
    for(int i=0; i<num_elements; i++){
    	if(fabs(data[i] - dec_data[i]) > max_err){
    		max_err = fabs(data[i] - dec_data[i]);
    	}
    }
    cout << "Max error = " << max_err << endl;
    free(dec_data);
}

template <class T>
void test(string filename){
    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    T max_val = 0;
    for(int i=0; i<num_elements; i++){
    	if(fabs(data[i]) > max_val) max_val = fabs(data[i]);
    }
    int level_exp = 0;
    frexp(max_val, &level_exp);

    evaluate<T, uint8_t>(data, num_elements, level_exp, 32);
    evaluate<T, uint16_t>(data, num_elements, level_exp, 32);
    evaluate<T, uint32_t>(data, num_elements, level_exp, 32);
    evaluate<T, uint64_t>(data, num_elements, level_exp, 32);
}

int main(int argc, char ** argv){

	string filename = string(argv[1]);
	test<float>(filename);
	return 0;

}