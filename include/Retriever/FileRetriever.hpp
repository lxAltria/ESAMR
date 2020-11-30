#ifndef _MDR_FILE_RETRIEVER_HPP
#define _MDR_FILE_RETRIEVER_HPP

#include "RetrieverInterface.hpp"
#include <cstdio>

namespace MDR {
    // Data retriever for files
    class FileRetriever : public concepts::RetrieverInterface {
    public:
        FileRetriever(const std::string& metadata_file, const std::vector<std::string>& level_files) : metadata_file(metadata_file), level_files(level_files) {}

        std::vector<uint8_t*> retrieve_level_components(const std::vector<uint32_t>& offsets, const std::vector<uint32_t>& retrieve_sizes) const {
            assert(offsets.size() == retrieve_sizes.size());
            std::vector<uint8_t*> level_components;
            for(int i=0; i<level_files.size(); i++){
                FILE * file = fopen(level_files[i].c_str(), "r");
                if(fseek(file, offsets[i], SEEK_SET)){
                    std::cerr << "Errors in fseek while retrieving from file" << std::endl;
                }
                uint8_t * buffer = (uint8_t *) malloc(retrieve_sizes[i]);
                fread(buffer, sizeof(uint8_t), retrieve_sizes[i], file);
                level_components.push_back(buffer);
            }
            return level_components;
        }

        uint8_t * load_metadata() const {
            FILE * file = fopen(metadata_file.c_str(), "r");
            fseek(file, 0, SEEK_END);
            uint32_t num_bytes = ftell(file);
            rewind(file);
            uint8_t * metadata = (uint8_t *) malloc(num_bytes);
            fread(metadata, 1, num_bytes, file);
            return metadata;
        }

        ~FileRetriever(){}

        void print() const {
            std::cout << "File retriever." << std::endl;
        }
    private:
        std::vector<std::string> level_files;
        std::string metadata_file;
    };
}
#endif
