#ifndef _LOSSLESS_UTILS_HPP
#define _LOSSLESS_UTILS_HPP

#include <vector>
#include <cstdlib>
#include <cstring>
#include "zstd.h"

namespace LOSSLESS{

using namespace std;

#define ZSTD_COMPRESSOR 1
uint32_t zstd_lossless_compress(int losslessCompressor, int level, uint8_t* data, uint32_t dataLength, uint8_t** compressBytes)
{
    uint32_t outSize = 0; 
    size_t estimatedCompressedSize = 0;
    switch(losslessCompressor)
    {
    case ZSTD_COMPRESSOR:
        if(dataLength < 100) 
            estimatedCompressedSize = 200;
        else
            estimatedCompressedSize = dataLength*1.2;
        *compressBytes = (uint8_t*)malloc(estimatedCompressedSize);
        *reinterpret_cast<size_t*>(*compressBytes) = dataLength;
        outSize = ZSTD_compress(*compressBytes + sizeof(size_t), estimatedCompressedSize, data, dataLength, level); //default setting of level is 3
        break;
    default:
        printf("Error: Unrecognized lossless compressor in sz_lossless_compress()\n");
    }
    return outSize + sizeof(size_t);
}
uint32_t zstd_lossless_decompress(int losslessCompressor, const uint8_t* compressBytes, uint32_t cmpSize, uint8_t** oriData)
{
    uint32_t outSize = 0;
    switch(losslessCompressor)
    {
    case ZSTD_COMPRESSOR:
        outSize = *reinterpret_cast<const size_t*>(compressBytes);
        *oriData = (uint8_t*)malloc(outSize);
        ZSTD_decompress(*oriData, outSize, compressBytes + sizeof(size_t), cmpSize - sizeof(size_t));
        break;
    default:
        printf("Error: Unrecognized lossless compressor in sz_lossless_decompress()\n");
    }
    return outSize;
}

}
#endif