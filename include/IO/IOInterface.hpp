#ifndef _MDR_IO_INTERFACE_HPP
#define _MDR_IO_INTERFACE_HPP

#include <cassert>

namespace MDR {
    namespace concepts {

        class IOInterface {
        public:

            virtual ~IOInterface() = default;

            virtual std::vector<std::vector<const uint8_t*>> retrieve_level_components(const std::vector<std::vector<uint32_t>>& level_sizes, const std::vector<uint32_t>& retrieve_sizes, const std::vector<uint8_t>& prev_level_num_bitplanes, const std::vector<uint8_t>& level_num_bitplanes) = 0;

            virtual uint8_t * load_metadata() const = 0;

            virtual void release() = 0;

            virtual void print() const = 0;

            virtual std::vector<std::vector<uint32_t>> write_level_components(const std::vector<std::vector<uint8_t*>>& level_components, const std::vector<std::vector<uint32_t>>& level_sizes) const = 0;

            virtual void write_metadata(uint8_t const * metadata, uint32_t size) const = 0;

            virtual void print() const = 0;

        };
    }
}
#endif
