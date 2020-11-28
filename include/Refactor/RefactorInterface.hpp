#ifndef _MDR_REFACTOR_INTERFACE_HPP
#define _MDR_REFACTOR_INTERFACE_HPP

namespace MDR {
    namespace concepts {

        // refactor: a general interface for scnetific data refactor
        template<class T>
        class RefactorInterface {
        public:

            virtual ~RefactorInterface() = default;

            virtual uint8_t * refactor(T const * data, const std::vector<uint32_t>& dims, uint8_t target_level, uint8_t num_bitplanes, uint32_t& refactored_size) = 0;

            virtual uint8_t * dump_metadata(uint32_t& metadata_size) const = 0;

            virtual void print() const = 0;
        };
    }
}
#endif
