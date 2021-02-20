#ifndef _MDR_ADAPTIVE_LEVEL_COMPRESSOR_HPP
#define _MDR_ADAPTIVE_LEVEL_COMPRESSOR_HPP

#include "LevelCompressorInterface.hpp"
#include "LosslessCompressor.hpp"

namespace MDR {
    #define LATTER_INDEX 24
    #define CR_THRESHOLD 1.1
    // compress all layers
    class AdaptiveLevelCompressor : public concepts::LevelCompressorInterface {
    public:
        AdaptiveLevelCompressor(){}
        void compress_level(std::vector<uint8_t*>& streams, std::vector<uint32_t>& stream_sizes) const {
            int stopping_index = 0;
            for(int i=0; i<streams.size(); i++){
                uint8_t * compressed = NULL;
                auto compressed_size = ZSTD::compress(streams[i], stream_sizes[i], &compressed);
                free(streams[i]);
                // skip the first
                if(i && (compressed_size / stream_sizes[i] < CR_THRESHOLD)){
                    stopping_index = i;
                    break;
                }
                streams[i] = compressed;
                stream_sizes[i] = compressed_size;
            }
            int latter_start_index = (stopping_index < LATTER_INDEX) ? LATTER_INDEX : stopping_index;
            for(int i=latter_start_index; i<streams.size(); i++){
                uint8_t * compressed = NULL;
                auto compressed_size = ZSTD::compress(streams[i], stream_sizes[i], &compressed);
                free(streams[i]);
                streams[i] = compressed;
                stream_sizes[i] = compressed_size;
            }
            std::cout << "stopping_index = " << stopping_index << std::endl;
        }
        void decompress_level(std::vector<const uint8_t*>& streams, const std::vector<uint32_t>& stream_sizes, uint8_t starting_bitplane, uint8_t num_bitplanes) {
            int stopping_index = 10;
            for(int i=0; i<num_bitplanes; i++){
                int bitplane_index = starting_bitplane + i;
                if((bitplane_index <= stopping_index) || (bitplane_index >= LATTER_INDEX)){
                    uint8_t * decompressed = NULL;
                    auto decompressed_size = ZSTD::decompress(streams[i], stream_sizes[bitplane_index], &decompressed);
                    buffer.push_back(decompressed);
                    streams[i] = decompressed;                    
                }
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
        ~AdaptiveLevelCompressor(){
            decompress_release();
        }
    private:
        std::vector<uint8_t*> buffer;
        std::vector<uint8_t> stopping_indices;
    };
}
#endif
