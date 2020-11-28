#ifndef _MDR_RETRIEVER_INTERFACE_HPP
#define _MDR_RETRIEVER_INTERFACE_HPP

namespace MDR {
    namespace concepts {

        // Error-controlled data retrieval
        class RetrieverInterface {
        public:

            virtual ~RetrieverInterface() = default;

            virtual uint32_t interpret_size(const std::vector<std::vector<uint32_t>>& level_sizes, const std::vector<std::vector<double>>& level_errors, const std::vector<uint8_t>& order, double tolerance, std::vector<uint8_t>& index) const = 0;

            virtual void print() const = 0;
        };
    }
}
#endif
