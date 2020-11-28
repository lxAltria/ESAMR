#ifndef _MDR_RECONSTRUCTOR_INTERFACE_HPP
#define _MDR_RECONSTRUCTOR_INTERFACE_HPP

namespace MDR {
    namespace concepts {

        // reconstructor: a general interface for scientific data reconstructor
        template<class T>
        class ReconstructorInterface {
        public:

            virtual ~ReconstructorInterface() = default;

            virtual T * reconstruct(uint8_t const * refactored_data, double tolerance) = 0;

            virtual void load_metadata(uint8_t const * metadata) = 0;

            virtual void print() const = 0;
        };
    }
}
#endif
