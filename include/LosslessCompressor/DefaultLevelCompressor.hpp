#ifndef _MDR_DEFAULT_LEVEL_COMPRESSOR_HPP
#define _MDR_DEFAULT_LEVEL_COMPRESSOR_HPP

#include "LevelCompressorInterface.hpp"
#include "LosslessCompressor.hpp"

namespace MDR {
    // compress all layers
    class DefaultLevelCompressor : public concepts::LevelCompressorInterface {
    public:
        DefaultLevelCompressor(){}
        void compress_level(std::vector<uint8_t*>& streams, std::vector<uint32_t>& stream_sizes) const {
            for(int i=0; i<streams.size(); i++){
                uint8_t * compressed = NULL;
                auto compressed_size = ZSTD::compress(streams[i], stream_sizes[i], &compressed);
                free(streams[i]);
                streams[i] = compressed;
                stream_sizes[i] = compressed_size;
            }
        }
        void decompress_level(std::vector<const uint8_t*>& streams, const std::vector<uint32_t>& stream_sizes, uint8_t starting_bitplane, uint8_t num_bitplanes) {
            for(int i=0; i<num_bitplanes; i++){
                uint8_t * decompressed = NULL;
                auto decompressed_size = ZSTD::decompress(streams[i], stream_sizes[starting_bitplane + i], &decompressed);
                buffer.push_back(decompressed);
                streams[i] = decompressed;
            }
        }
        void decompress_release(){
            for(int i=0; i<buffer.size(); i++){
                free(buffer[i]);
            }
            buffer.clear();
        }
        void print() const {
            std::cout << "Default Level lossless compressor" << std::endl;
        }
        ~DefaultLevelCompressor(){
            decompress_release();
        }
    private:
        std::vector<uint8_t*> buffer;
    };
}
#endif
