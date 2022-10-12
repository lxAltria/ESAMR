#ifndef _MDR_FILE_WRITER_HPP
#define _MDR_FILE_WRITER_HPP

#include "WriterInterface.hpp"
#include <cstdio>

namespace MDR {

    // A writer that writes the concatenated level components
    class ConcatLevelFileWriter : public concepts::WriterInterface {
    public:
        ConcatLevelFileWriter(const std::string& metadata_file, const std::vector<std::string>& level_files) : metadata_file(metadata_file), level_files(level_files) {}

        std::vector<uint32_t> write_level_components(const std::vector<std::vector<uint8_t*>>& level_components, const std::vector<std::vector<uint32_t>>& level_sizes) const {
            std::vector<uint32_t> level_num;
            for(int i=0; i<level_components.size(); i++){
                uint32_t concated_level_size = 0;
                for(int j=0; j<level_components[i].size(); j++){
                    concated_level_size += level_sizes[i][j];
                }
                uint8_t * concated_level_data = (uint8_t *) malloc(concated_level_size);
                uint8_t * concated_level_data_pos = concated_level_data;
                for(int j=0; j<level_components[i].size(); j++){
                    memcpy(concated_level_data_pos, level_components[i][j], level_sizes[i][j]);
                    concated_level_data_pos += level_sizes[i][j];
                }
                FILE * file = fopen((level_files[i]).c_str(), "w");
                fwrite(concated_level_data, 1, concated_level_size, file);
                fclose(file);
                free(concated_level_data);
                level_num.push_back(1);
            }
            return level_num;
        }

        void write_metadata(uint8_t const * metadata, uint32_t size) const {
            FILE * file = fopen(metadata_file.c_str(), "w");
            fwrite(metadata, 1, size, file);
            fclose(file);
        }

        ~ConcatLevelFileWriter(){}

        void print() const {
            std::cout << "File writer." << std::endl;
        }
    private:
        std::vector<std::string> level_files;
        std::string metadata_file;
    };

    // A writer that writes individual level components
    class DirectFileWriter : public concepts::WriterInterface {
    public:
        DirectFileWriter(const std::string& metadata_file, const std::string& data_name) : metadata_file(metadata_file), data_name(data_name) {}

        std::vector<uint32_t> write_level_components(const std::vector<std::vector<uint8_t*>>& level_components, const std::vector<std::vector<uint32_t>>& level_sizes) const {
            std::vector<uint32_t> level_num;
            for(int i=0; i<level_components.size(); i++){
                level_num.push_back(level_components[i].size());
                for(int j=0; j<level_components[i].size(); j++){
                    std::string filename = data_name + "_" + std::to_string(i) + "_" + std::to_string(j);
                    FILE * file = fopen(filename.c_str(), "w");
                    fwrite(level_components[i][j], 1, level_sizes[i][j], file);
                    fclose(file);
                }
            }
            return level_num;
        }

        void write_metadata(uint8_t const * metadata, uint32_t size) const {
            FILE * file = fopen(metadata_file.c_str(), "w");
            fwrite(metadata, 1, size, file);
            fclose(file);
        }

        ~DirectFileWriter(){}

        void print() const {
            std::cout << "Direct file writer." << std::endl;
        }
    protected:
        std::string metadata_file;
        std::string data_name;
    };

    // A writer that writes multiblock components
    // used in adaptive refactor
    class MultiblockFileWriter : public DirectFileWriter {
    public:
        MultiblockFileWriter(const std::string& metadata_file, const std::string& data_name, int num_bitplanes=32) : DirectFileWriter(metadata_file, data_name), num_bitplanes(num_bitplanes) {}

        std::vector<uint32_t> write_level_components(const std::vector<std::vector<uint8_t*>>& level_components, const std::vector<std::vector<uint32_t>>& level_sizes, const std::vector<std::vector<int>>& level_merge_counts) const {
            std::vector<uint32_t> level_num;
            for(int i=0; i<level_components.size(); i++){
                level_num.push_back(level_components[i].size());
                int bp_num = 0;
                int count = 0;
                int block_num = 0;
                for(int j=0; j<level_components[i].size(); j++){
                    std::string filename = data_name + "_L" + std::to_string(i) + "_B" + std::to_string(block_num) + "_S" + std::to_string(count);
                    FILE * file = fopen(filename.c_str(), "w");
                    fwrite(level_components[i][j], 1, level_sizes[i][j], file);
                    fclose(file);
                    count ++;
                    bp_num += level_merge_counts[i][j];
                    if(bp_num % num_bitplanes == 0){
                        block_num ++;
                        count = 0;
                    }
                }
            }
            return level_num;
        }

        ~MultiblockFileWriter(){}

        void print() const {
            std::cout << "Multiblock file writer." << std::endl;
        }
    private:
        int num_bitplanes;
    };
}
#endif
