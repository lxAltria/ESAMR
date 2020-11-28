#ifndef _MDR_LOSSLESS_COMPRESSOR_INTERFACE_HPP
#define _MDR_LOSSLESS_COMPRESSOR_INTERFACE_HPP

namespace MDR {
    namespace concepts {

        // interface for lossless compressor 
        class LosslessCompressorInterface {
        public:

            virtual ~LosslessCompressorInterface() = default;

            virtual uint32_t compress(uint8_t* data, uint32_t dataLength, uint8_t** compressBytes) const  = 0;

            virtual uint32_t decompress(const uint8_t* compressBytes, uint32_t cmpSize, uint8_t** oriData) const = 0;

            virtual void print() const = 0;
        };
    }
}
#endif
