#ifndef _MDR_NULL_LOSSLESS_COMPRESSOR_HPP
#define _MDR_NULL_LOSSLESS_COMPRESSOR_HPP

#include "LosslessCompressorInterface.hpp"

namespace MDR {
    // Null lossless compressor
    class NullLosslessCompressor : public concepts::LosslessCompressorInterface {
    public:
        NullLosslessCompressor(){}
        uint32_t compress(uint8_t* data, uint32_t dataLength, uint8_t** compressBytes) const {
            *compressBytes = data;
            return dataLength;
        }
        uint32_t decompress(const uint8_t* compressBytes, uint32_t cmpSize, uint8_t** oriData) const {
            *oriData = compressBytes;
            return cmpSize;
        }
        void print() const {
            std::cout << "Null lossless compressor" << std::endl;
        }
    };
}
#endif
