#include<iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include "utils.hpp"
#include "encode.hpp"
#include "decode.hpp"

using namespace std;

// encode n < 32 integer data
// template <class T>
uint32_t encode_block(uint32_t * data, int n, vector<uint32_t *>& encoded, vector<uint32_t>& offset){
	bool recorded = false;
	uint32_t block_indicator = 0;
	int prec = 31;
	for(int k=prec; k>=0; k--){
		uint32_t bitplane_value = 0;
		uint32_t bitplane_index = prec - k;
		for (int i=0; i<n; i++){
			bitplane_value += (uint32_t)((data[i] >> k) & 1u) << i;
		}
		if(bitplane_value || recorded){
			if(!recorded){
				block_indicator = bitplane_index;
				recorded = true;
			}
			encoded[bitplane_index][offset[bitplane_index] ++] = bitplane_value;
		}
	}
	return block_indicator;
}

template <class T>
vector<uint32_t *> encode_data(T * data, int n, int level_exp){
	uint8_t block_size = 32;
	uint32_t int_data[32] = {0};
	const int prec = 32;
	vector<uint32_t *> encode_bitplanes;
	vector<uint32_t> offset = vector<uint32_t>(prec, 0);
	for(int i=0; i<prec; i++){
		encode_bitplanes.push_back((uint32_t *) malloc(n + sizeof(uint32_t)));
	}
	for(int i=0; i<n; i+=block_size){
		for(int j=0; j<block_size; j++){
			T cur_data = ldexp(data[i + j], prec - 1 - level_exp);
			int64_t fix_point = (int64_t) cur_data;
			bool sign = data[i + j] < 0;
			int_data[j] = sign ? -fix_point : +fix_point;
		}
		encode_block(int_data, block_size, encode_bitplanes, offset);
	}
	int sum = 0;
	for(int i=0; i<offset.size(); i++){
		sum += offset[i];
		cout << offset[i] << " ";
	}
	cout << endl;
	cout << sum << endl;
	return encode_bitplanes;
}

template <class T, class T2>
void evaluate(const vector<T>& data, size_t num_elements, int level_exp, int num_bitplanes){

    struct timespec start, end;
    int err = 0;
    err = clock_gettime(CLOCK_REALTIME, &start);
    vector<uint8_t> sbp;
    vector<uint32_t> sizes;
    auto result = REFACTOR::encode<T, T2>(data.data(), num_elements, level_exp, num_bitplanes, sbp, sizes);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << result.size() << endl;
    cout << "Encoding time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;	

    // err = clock_gettime(CLOCK_REALTIME, &start);
    // auto dec_data = REFACTOR::decode<T, T2>(result, sbp, num_elements, level_exp, num_bitplanes);
    // err = clock_gettime(CLOCK_REALTIME, &end);
    // cout << "Decoding time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;

    // cout << "Encoded sizes: ";
    // for(int i=0; i<sizes.size(); i++){
    // 	cout << sizes[i] << " ";
    // }
    // cout << endl;

    // T max_err = 0;
    // for(int i=0; i<num_elements; i++){
    // 	if(fabs(data[i] - dec_data[i]) > max_err){
    // 		max_err = fabs(data[i] - dec_data[i]);
    // 		if(max_err > 113) cout << i << ": " << data[i] << " " << dec_data[i] << endl;
    // 	}
    // }
    // cout << "Max error = " << max_err << endl;
    // free(dec_data);
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
    // struct timespec start, end;
    // int err = 0;
    // err = clock_gettime(CLOCK_REALTIME, &start);
    // auto res = encode_data(data.data(), num_elements, level_exp);
    // err = clock_gettime(CLOCK_REALTIME, &end);
    // cout << res.size() << endl;
    // cout << "Encoding time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    // err = clock_gettime(CLOCK_REALTIME, &start);

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