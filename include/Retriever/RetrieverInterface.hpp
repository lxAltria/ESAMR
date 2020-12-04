#ifndef _MDR_RETRIEVER_INTERFACE_HPP
#define _MDR_RETRIEVER_INTERFACE_HPP

#include <cassert>

namespace MDR {
    namespace concepts {

        // Error-controlled data retrieval
        class RetrieverInterface {
        public:

            virtual ~RetrieverInterface() = default;

            virtual std::vector<uint8_t*> retrieve_level_components(const std::vector<uint32_t>& offsets, const std::vector<uint32_t>& retrieve_sizes) const = 0;

            virtual uint8_t * load_metadata() const = 0;

            virtual void print() const = 0;
        };
    }
}
#endif
