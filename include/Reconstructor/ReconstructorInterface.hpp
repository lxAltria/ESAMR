#ifndef _MDR_RECONSTRUCTOR_INTERFACE_HPP
#define _MDR_RECONSTRUCTOR_INTERFACE_HPP

namespace MDR {
    namespace concepts {

        // reconstructor: a general interface for scientific data reconstructor
        template<class T>
        class ReconstructorInterface {
        public:

            virtual ~ReconstructorInterface() = default;

            virtual T * reconstruct(uint8_t const * retrieved_data, uint32_t retrieved_size) = 0;

            virtual T * progressive_reconstruct(uint8_t const * retrieved_data, uint32_t retrieved_size) = 0;

            virtual void load_metadata(uint8_t const * metadata) = 0;

            virtual void print() const = 0;
        };
    }
}
#endif
