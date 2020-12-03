#ifndef _MDR_PERBIT_BP_ENCODER_HPP
#define _MDR_PERBIT_BP_ENCODER_HPP

#include "BitplaneEncoderInterface.hpp"

namespace MDR {
    class BitEncoder{
    public:
        BitEncoder(uint64_t * stream_begin_pos){
            stream_begin = stream_begin_pos;
            stream_pos = stream_begin;
            buffer = 0;
            position = 0;
        }
        void encode(uint8_t b){
            buffer += b << (position++);
            if(position == 64){
                *(stream_pos ++) = buffer;
                buffer = 0;
                position = 0;
            }
        }
        void flush(){
            if(position){
                *(stream_pos ++) = buffer;
                buffer = 0;
                position = 0;
            }
        }
        uint32_t size(){
            return (stream_pos - stream_begin);
        }
    private:
        uint64_t buffer = 0;
        uint8_t position = 0;
        uint64_t * stream_pos = NULL;
        uint64_t * stream_begin = NULL;
    };

    class BitDecoder{
    public:
        BitDecoder(uint64_t * stream_begin_pos){
            stream_begin = stream_begin_pos;
            stream_pos = stream_begin;
            buffer = *(stream_pos ++);
            position = 64;
        }
        uint8_t decode(){
            uint8_t b = buffer & 1u;
            buffer >>= 1;
            position --;
            if(position == 0){
                buffer = *(stream_pos ++);
                position = 64;
            }
            return b;
        }
    private:
        uint64_t buffer = 0;
        uint8_t position = 0;
        uint64_t const * stream_pos = NULL;
        uint64_t const * stream_begin = NULL;
    };


    // per bit bitplane encoder that encodes data by bit using T_stream type buffer
    template<class T_data, class T_stream>
    class PerBitBPEncoder : public concepts::BitplaneEncoderInterface<T_data> {
    public:
        PerBitBPEncoder(){
            static_assert(std::is_floating_point<T_data>::value, "PerBitBPEncoder: input data must be floating points.");
            static_assert(!std::is_same<T_data, long double>::value, "PerBitBPEncoder: long double is not supported.");
            static_assert(std::is_unsigned<T_stream>::value, "PerBitBPEncoder: streams must be unsigned integers.");
            static_assert(std::is_integral<T_stream>::value, "PerBitBPEncoder: streams must be unsigned integers.");
        }

        std::vector<uint8_t *> encode(T_data const * data, int32_t n, int32_t exp, uint8_t num_bitplanes, std::vector<uint32_t>& stream_sizes) const {
            assert(num_bitplanes > 0);
            // determine block size based on bitplane integer type
            const uint32_t block_size = 64;
            std::vector<uint8_t> starting_bitplanes = std::vector<uint8_t>((n - 1)/block_size + 1, 0);
            stream_sizes = std::vector<uint32_t>(num_bitplanes, 0);
            // define fixed point type
            using T_fp = typename std::conditional<std::is_same<T_data, double>::value, uint64_t, uint32_t>::type;
            std::vector<uint8_t *> streams;
            for(int i=0; i<num_bitplanes; i++){
                streams.push_back((uint8_t *) malloc(2 * n / UINT8_BITS + sizeof(T_stream)));
            }
            std::vector<T_fp> int_data_buffer(block_size, 0);
            std::vector<uint8_t> signs(block_size, 0);
            std::vector<T_stream *> streams_pos(streams.size());
            for(int i=0; i<streams.size(); i++){
                streams_pos[i] = reinterpret_cast<T_stream*>(streams[i]);
            }
            T_data const * data_pos = data;
            for(int i=0; i<n - block_size; i+=block_size){
                T_stream sign_bitplane = 0;
                for(int j=0; j<block_size; j++){
                    T_data cur_data = *(data_pos++);
                    T_data shifted_data = ldexp(cur_data, num_bitplanes - exp);
                    int64_t fix_point = (int64_t) shifted_data;
                    bool sign = cur_data < 0;
                    int_data_buffer[j] = sign ? -fix_point : +fix_point;
                    signs[j] = sign;
                }
                encode_block(int_data_buffer.data(), signs.data(), block_size, num_bitplanes, streams_pos);
            }
            // leftover
            {
                int rest_size = n - block_size * block_id;
                T_stream sign_bitplane = 0;
                for(int j=0; j<rest_size; j++){
                    T_data cur_data = *(data_pos++);
                    T_data shifted_data = ldexp(cur_data, num_bitplanes - exp);
                    int64_t fix_point = (int64_t) shifted_data;
                    T_stream sign = cur_data < 0;
                    int_data_buffer[j] = sign ? -fix_point : +fix_point;
                    sign_bitplane += sign << j;
                }
                starting_bitplanes[block_id ++] = encode_block(int_data_buffer.data(), rest_size, num_bitplanes, sign_bitplane, streams_pos);
            }
            for(int i=0; i<num_bitplanes; i++){
                stream_sizes[i] = reinterpret_cast<uint8_t*>(streams_pos[i]) - streams[i];
            }
            return streams;
        }

        T_data * decode(const std::vector<uint8_t const *>& streams, size_t n, int exp, uint8_t num_bitplanes) const {
            assert(num_bitplanes > 0);
            uint32_t block_size = block_size_based_on_bitplane_int_type<T_stream>();
            // define fixed point type
            using T_fp = typename std::conditional<std::is_same<T_data, double>::value, uint64_t, uint32_t>::type;
            // define manttisa
            T_data * data = (T_data *) malloc(n * sizeof(T_data));
            std::vector<T_stream const *> streams_pos(streams.size());
            for(int i=0; i<streams.size(); i++){
                streams_pos[i] = reinterpret_cast<T_stream const *>(streams[i]);
            }
            // deinterleave the first bitplane
            uint32_t starting_bitplane_size = *reinterpret_cast<int32_t const*>(streams_pos[0]);
            uint8_t const * starting_bitplanes = reinterpret_cast<uint8_t const*>(streams_pos[0]) + sizeof(uint32_t);
            streams_pos[0] = reinterpret_cast<T_stream const *>(starting_bitplanes + starting_bitplane_size);

            std::vector<T_fp> int_data_buffer(block_size, 0);
            // decode
            T_data * data_pos = data;
            int block_id = 0;
            for(int i=0; i<n - block_size; i+=block_size){
                memset(int_data_buffer.data(), 0, block_size * sizeof(T_fp));
                uint8_t starting_bitplane = starting_bitplanes[block_id ++];
                T_stream sign_bitplane = 0;
                if(starting_bitplane < num_bitplanes){
                    sign_bitplane = decode_block(streams_pos, block_size, num_bitplanes - starting_bitplane, starting_bitplane, int_data_buffer.data());
                }
                for(int j=0; j<block_size; j++, sign_bitplane >>= 1){
                    T_data cur_data = ldexp((T_data)int_data_buffer[j], - num_bitplanes + exp);
                    *(data_pos++) = (sign_bitplane & 1u) ? -cur_data : cur_data;
                }
            }
            // leftover
            {
                int rest_size = n - block_size * block_id;
                memset(int_data_buffer.data(), 0, block_size * sizeof(T_fp));
                int starting_bitplane = starting_bitplanes[block_id ++];
                T_stream sign_bitplane = 0;
                if(starting_bitplane < num_bitplanes){
                    sign_bitplane = decode_block(streams_pos, rest_size, num_bitplanes - starting_bitplane, starting_bitplane, int_data_buffer.data());
                }
                for(int j=0; j<rest_size; j++, sign_bitplane >>= 1){
                    T_data cur_data = ldexp((T_data)int_data_buffer[j], - num_bitplanes + exp);
                    *(data_pos++) = (sign_bitplane & 1u) ? -cur_data : cur_data;
                }
            }
            return data;
        }

        void print() const {
            std::cout << "Per-bit bitplane encoder" << std::endl;
        }
    private:
        template <class T_int>
        uint8_t encode_block(T_int const * data, size_t n, uint8_t num_bitplanes, T_stream sign, std::vector<T_stream *>& streams_pos) const {
            bool recorded = false;
            uint8_t starting_bitplane = num_bitplanes;
            for(int k=num_bitplanes - 1; k>=0; k--){
                T_stream bitplane_value = 0;
                T_stream bitplane_index = num_bitplanes - 1 - k;
                for (int i=0; i<n; i++){
                    bitplane_value += (T_stream)((data[i] >> k) & 1u) << i;
                }
                if(bitplane_value || recorded){
                    if(!recorded){
                        recorded = true;
                        starting_bitplane = bitplane_index;
                        *(streams_pos[bitplane_index] ++) = sign;
                    }
                    *(streams_pos[bitplane_index] ++) = bitplane_value;
                }
            }
            return starting_bitplane;
        }

        template <class T_int>
        T_stream decode_block(std::vector<T_stream const *>& streams_pos, size_t n, uint8_t num_bitplanes, uint8_t starting_bitplane, T_int * data) const {
            T_stream sign_bitplane = *(streams_pos[starting_bitplane] ++);
            for(int k=num_bitplanes - 1; k>=0; k--){
                T_stream bitplane_index = starting_bitplane + num_bitplanes - 1 - k;
                T_stream bitplane_value = *(streams_pos[bitplane_index] ++);
                for (int i=0; i<n; i++){
                    data[i] += ((bitplane_value >> i) & 1u) << k;
                }
            }
            return sign_bitplane;
        }

    };
}
#endif
