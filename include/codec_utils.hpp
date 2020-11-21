#ifndef _REFACTOR_CODEC_UTILS_HPP
#define _REFACTOR_CODEC_UTILS_HPP

#include <iostream>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <bitset>
#include <iomanip>
#include <climits>
#include <type_traits>

namespace REFACTOR{

// number of bits in a char
#define UINT8_BITS 8


template<class T>
inline int block_size_based_on_bitplane_int_type(){
	int block_size = 0;
	if(std::is_same<T, uint64_t>::value){
		block_size = 64;
	}
	else if(std::is_same<T, uint32_t>::value){
		block_size = 32;
	}
	else if(std::is_same<T, uint16_t>::value){
		block_size = 16;
	}
	else if(std::is_same<T, uint8_t>::value){
		block_size = 8;
	}
	else{
		std::cerr << "Integer type not supported." << std::endl;
		exit(0);
	}
	return block_size;
}

inline uint32_t compute_compact_size(uint32_t size, uint8_t index_size){
	return (size * index_size - 1) / UINT8_BITS + 1;
}

/* 
@params encoded: addresses of encoded bitplanes
@params offset: position of encoded bitplanes
return compact array
*/
template<class T>
uint8_t * compact(const std::vector<T>& data, uint8_t index_size){
	static_assert(std::is_unsigned<T>::value, "codec_utils compact: input array must be unsigned integers.");
	static_assert(std::is_integral<T>::value, "codec_utils compact: input array must be unsigned integers.");
	const uint8_t mask[UINT8_BITS + 1] = {0x00, 0x01, 0x03, 0x07, 0x0f, 0x1f, 0x3f, 0x7f, 0xff};
	uint32_t compact_size = compute_compact_size(data.size(), index_size);
	uint8_t * compact_data = (uint8_t *) malloc(compact_size);
	uint8_t * compact_data_pos = compact_data;
	uint8_t buffer = 0;
	uint8_t rest_bits = UINT8_BITS;
	for(int i=0; i<data.size(); i++){
		if(index_size <= rest_bits){
			rest_bits -= index_size;
			buffer <<= rest_bits;
			buffer += data[i];
		}
		else{
			buffer <<= rest_bits;
			buffer += data[i] >> (index_size - rest_bits);
			*(compact_data_pos ++) = buffer;
			buffer = data[i] & mask[index_size - rest_bits];
			rest_bits = UINT8_BITS + rest_bits - index_size;
		}
	}
	std::cout << std::endl;
	// flush buffer
	if(buffer){
		*(compact_data_pos ++) = buffer << rest_bits;
	}
	return compact_data;
}

template<class T>
std::vector<T> decompact(const uint8_t * compact_data, uint32_t size, uint8_t index_size){
	static_assert(std::is_unsigned<T>::value, "codec_utils decompact: input array must be unsigned integers.");
	static_assert(std::is_integral<T>::value, "codec_utils decompact: input array must be unsigned integers.");
	std::vector<T> data(size, 0);
	const uint8_t mask[UINT8_BITS + 1] = {0x00, 0x01, 0x03, 0x07, 0x0f, 0x1f, 0x3f, 0x7f, 0xff};
	const uint8_t * compact_data_pos = compact_data;
	uint8_t buffer = 0;
	uint8_t rest_bits = 0;
	for(int i=0; i<size; i++){
		if(index_size <= rest_bits){
			rest_bits -= index_size;
			data[i] = buffer >> rest_bits;
			buffer = buffer & mask[rest_bits];
		}
		else{
			data[i] = buffer << (index_size - rest_bits);
			buffer = *(compact_data_pos ++);
			rest_bits = UINT8_BITS + rest_bits - index_size;
			data[i] += buffer >> rest_bits;
			buffer = buffer & mask[rest_bits];
		}
	}
	return data;
}

}
#endif