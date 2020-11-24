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
#define UINT32_BITS 32
#define UINT64_BITS 64

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

}
#endif