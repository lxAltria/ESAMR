#ifndef _REFACTOR_ENCODE_HPP
#define _REFACTOR_ENCODE_HPP

#include "codec_utils.hpp"

using namespace std;

namespace REFACTOR{

/* 
T_fp: integer type for intermediate fixed point data, related to original data type
T_bitplane_int: integer type for encoded bitplane, related to block size
*/
// encode num_elements < 64 integer data
// uint64_t: num_elements < 64
// uint32_t: num_elements < 32
// uint16_t: num_elements < 16
// uint8_t: num_elements < 8
/*
@params int_data: int data to be encoded
@params num_elements: number of elements
@params starting_bitplane: the index of the first bitplane to be encoded
@params num_bitplanes: number of bitplanes to be encoded
@params encoded: addresses of encoded bitplanes
@params offset: position of encoded bitplanes
*/
template <class T_fp, class T_bitplane_int>
inline uint8_t encode_int_block_level_exp(const T_fp * const int_data, int num_elements, int num_bitplanes, T_bitplane_int sign, vector<T_bitplane_int *>& encoded_pos){
	bool recorded = false;
	uint8_t starting_bitplane = 64;
	for(int k=num_bitplanes - 1; k>=0; k--){
		T_bitplane_int bitplane_value = 0;
		T_bitplane_int bitplane_index = num_bitplanes - 1 - k;
		for (int i=0; i<num_elements; i++){
			bitplane_value += (T_bitplane_int)((int_data[i] >> k) & 1u) << i;
		}
		if(bitplane_value || recorded){
			if(!recorded){
				recorded = true;
				starting_bitplane = bitplane_index;
				*(encoded_pos[bitplane_index] ++) = sign;
			}
			*(encoded_pos[bitplane_index] ++) = bitplane_value;
		}
	}
	return starting_bitplane;
}

template<class T, class T_bitplane_int>
vector<T_bitplane_int *> encode(T const * data, size_t n, int aligned_exp, int num_bitplanes, vector<uint8_t>& starting_bitplanes, vector<uint32_t>& encoded_sizes){
	// determine block size based on bitplane integer type
	int block_size = block_size_based_on_bitplane_int_type<T_bitplane_int>();
	cout << "block_size = " << block_size << endl;
	starting_bitplanes = vector<uint8_t>((n - 1)/block_size + 1, 0);
	encoded_sizes = vector<uint32_t>(num_bitplanes, 0);
	// define fixed point type
	using T_fp = typename std::conditional<std::is_same<T, double>::value, uint64_t, uint32_t>::type;
	// define manttisa
	const int mantissa_bits = std::is_same<T_bitplane_int, double>::value ? 52 : 23;
	// vector<uint32_t> offset(num_bitplanes, 0);
	vector<T_bitplane_int *> encoded_bitplanes;
	for(int i=0; i<num_bitplanes; i++){
		encoded_bitplanes.push_back((T_bitplane_int *) malloc(2 * n + sizeof(T_bitplane_int)));
	}
	vector<T_fp> int_data_buffer(block_size, 0);
	vector<T_bitplane_int *> encoded_pos(encoded_bitplanes);
	T const * data_pos = data;
	int block_id=0;
	for(int i=0; i<n - block_size; i+=block_size){
		T_bitplane_int sign_bitplane = 0;
		for(int j=0; j<block_size; j++){
			T cur_data = *(data_pos++);
			T shifted_data = ldexp(cur_data, num_bitplanes - aligned_exp);
			int64_t fix_point = (int64_t) shifted_data;
			T_bitplane_int sign = cur_data < 0;
			int_data_buffer[j] = sign ? -fix_point : +fix_point;
			sign_bitplane += sign << j;
		}
		starting_bitplanes[block_id ++] = encode_int_block_level_exp(int_data_buffer.data(), block_size, num_bitplanes, sign_bitplane, encoded_pos);
	}
	// leftover
	{
		int rest_size = n - block_size * block_id;
		cout << rest_size << endl;
		// memset(int_data_buffer.data(), 0, block_size * sizeof(T_fp));
		T_bitplane_int sign_bitplane = 0;
		for(int j=0; j<rest_size; j++){
			T cur_data = *(data_pos++);
			T shifted_data = ldexp(cur_data, num_bitplanes - aligned_exp);
			int64_t fix_point = (int64_t) shifted_data;
			T_bitplane_int sign = cur_data < 0;
			int_data_buffer[j] = sign ? -fix_point : +fix_point;
			sign_bitplane += sign << j;
		}
		starting_bitplanes[block_id ++] = encode_int_block_level_exp(int_data_buffer.data(), rest_size, num_bitplanes, sign_bitplane, encoded_pos);
	}
	for(int i=0; i<num_bitplanes; i++){
		encoded_sizes[i] = (encoded_pos[i] - encoded_bitplanes[i]) * sizeof(T_bitplane_int);
	}
	return encoded_bitplanes;
}

}
#endif