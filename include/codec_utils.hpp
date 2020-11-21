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
	return (size * UINT8_BITS - 1) / index_size + 1;
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
	uint32_t compact_size = data.size() * sizeof(T);//compute_compact_size(data.size(), index_size);
	uint8_t * compact_data = (uint8_t *) malloc(compact_size);
	memcpy(compact_data, data.data(), data.size() * sizeof(T));
	return compact_data;
}

template<class T>
std::vector<T> decompact(const uint8_t * compact_data, uint32_t size){
	std::vector<T> data(size);
	memcpy(data.data(), compact_data, size * sizeof(T));
	return data;
}

}
#endif