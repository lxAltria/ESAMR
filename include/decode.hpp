#ifndef _REFACTOR_DECODE_HPP
#define _REFACTOR_DECODE_HPP

#include "codec_utils.hpp"

using namespace std;

namespace REFACTOR{

// encode num_elements < 64 32-bit integer data
// uint64_t: num_elements < 64
// uint32_t: num_elements < 32
// uint16_t: num_elements < 16
// uint8_t: num_elements < 8
/*
@params encoded: addresses of encoded bitplanes
@params offset: position of encoded bitplanes
@params num_elements: number of elements
@params starting_bitplane: the index of the first bitplane to be encoded
@params num_bitplanes: number of bitplanes to be encoded
@params data: decoded int data
*/
template <class T, class T_fp>
inline void decode_int_block(vector<const T *>& encoded_pos, int num_elements, int starting_bitplane, int num_bitplanes, T_fp * data){
	for(int k=num_bitplanes - 1; k>=0; k--){
		T bitplane_index = starting_bitplane + num_bitplanes - 1 - k;
		T bitplane_value = *(encoded_pos[bitplane_index] ++);
		for (int i=0; i<num_elements; i++){
			data[i] += ((bitplane_value >> i) & 1u) << k;
		}
	}
}

template<class T, class T_bitplane_int>
T * decode(const vector<const T_bitplane_int *>& encoded, const vector<uint8_t>& starting_bitplanes, size_t n, int aligned_exp, int num_bitplanes){
	const int block_size = block_size_based_on_bitplane_int_type<T_bitplane_int>();
	cout << "block_size = " << block_size << endl;
	// define fixed point type
	using T_fp = typename std::conditional<std::is_same<T, double>::value, uint64_t, uint32_t>::type;
	// define manttisa
	const int mantissa_bits = std::is_same<T_bitplane_int, double>::value ? 52 : 23;
	T * data = (T *) malloc(n * sizeof(T));
	vector<const T_bitplane_int *> encoded_pos(encoded);
	vector<T_fp> int_data_buffer(block_size, 0);
	// decode
	T * data_pos = data;
	int block_id = 0;
	for(int i=0; i<n - block_size; i+=block_size){
		memset(int_data_buffer.data(), 0, block_size * sizeof(T_fp));
		int starting_bitplane = starting_bitplanes[block_id ++];
		T_bitplane_int sign_bitplane = 0;
		if(starting_bitplane < num_bitplanes){
	        sign_bitplane = *(encoded_pos[starting_bitplane] ++);
	        decode_int_block(encoded_pos, block_size, starting_bitplane, num_bitplanes - starting_bitplane, int_data_buffer.data());
    	}
        for(int j=0; j<block_size; j++, sign_bitplane >>= 1){
	        T cur_data = ldexp((T)int_data_buffer[j], - num_bitplanes + aligned_exp);
	        *(data_pos++) = (sign_bitplane & 1u) ? -cur_data : cur_data;
    	}
    }
    // leftover
    {
        int rest_size = n - block_size * block_id;
        memset(int_data_buffer.data(), 0, block_size * sizeof(T_fp));
        int starting_bitplane = starting_bitplanes[block_id ++];
        T_bitplane_int sign_bitplane = 0;
        if(starting_bitplane < num_bitplanes){
            sign_bitplane = *(encoded_pos[starting_bitplane] ++);
            decode_int_block(encoded_pos, rest_size, starting_bitplane, num_bitplanes - starting_bitplane, int_data_buffer.data());
        }
        for(int j=0; j<rest_size; j++, sign_bitplane >>= 1){
            T cur_data = ldexp((T)int_data_buffer[j], - num_bitplanes + aligned_exp);
            *(data_pos++) = (sign_bitplane & 1u) ? -cur_data : cur_data;
        }
    }
    return data;
}

}
#endif