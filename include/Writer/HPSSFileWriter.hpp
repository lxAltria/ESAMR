#ifndef _MDR_HPSS_WRITER_HPP
#define _MDR_HPSS_WRITER_HPP

#include "WriterInterface.hpp"
#include <cstdio>
#include <mpi.h>
#include <adios2.h>

namespace MDR {
    // A writer that writes the concatenated level components
    // Merge multiple components if size is small
    class HPSSFileWriter : public concepts::WriterInterface {
    public:
        HPSSFileWriter(const std::string& metadata_file, const std::vector<std::string>& level_files, int num_process, int min_HPSS_size) : metadata_file(metadata_file), level_files(level_files), min_size((min_HPSS_size - 1)/num_process + 1) {}

        std::vector<std::vector<uint32_t>> write_level_components(const std::vector<std::vector<uint8_t*>>& level_components, const std::vector<std::vector<uint32_t>>& level_sizes) const {
            int rank, size;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_Comm_size(MPI_COMM_WORLD, &size);
            adios2::ADIOS ad(MPI_COMM_WORLD);
            adios2::IO bpIO = ad.DeclareIO("WritePrecisionSegments");

            std::vector<std::vector<uint32_t>> level_merged_count;
            for(int i=0; i<level_components.size(); i++){
                std::vector<uint32_t> merged_count;
                uint64_t concated_level_size = 0;
                uint32_t prev_index = 0;
                uint32_t count = 0;
                for(int j=0; j<level_components[i].size(); j++){
                    concated_level_size += level_sizes[i][j];
                    // concate bitplane in [prev_index, j]
                    if((concated_level_size >= min_size) || (j == level_components[i].size() - 1)){
                        int num_bitplanes = j - prev_index + 1;
                        merged_count.push_back(num_bitplanes);
                        // TODO: deal with the last file that may not be larger than min_size
                        uint8_t * concated_level_data = (uint8_t *) malloc(concated_level_size);
                        uint8_t * concated_level_data_pos = concated_level_data;
                        // std::cout << +prev_index << " " << j << " " << concated_level_size << std::endl;
                        for(int k=prev_index; k<=j; k++){
                            memcpy(concated_level_data_pos, level_components[i][k], level_sizes[i][k]);
                            concated_level_data_pos += level_sizes[i][k];
                        }
                        // FILE * file = fopen((level_files[i] + "_" + std::to_string(count)).c_str(), "w");
                        // fwrite(concated_level_data, 1, concated_level_size, file);
                        // fclose(file);
                        // MPI_IO
                        std::string filename = level_files[i] + "_" + std::to_string(count);
                        // MPIIO_write(MPI_COMM_WORLD, rank, size, concated_level_size, concated_level_data, filename);
                        {
                            adios2::Variable<uint8_t> bp_fdata = bpIO.DefineVariable<uint8_t>(
                                  filename, {(uint64_t)concated_level_size * size}, {(uint64_t)concated_level_size * rank}, {(uint64_t)concated_level_size}, adios2::ConstantDims);
                            adios2::Engine bpFileWriter = bpIO.Open(filename, adios2::Mode::Write);
                            bpFileWriter.Put<uint8_t>(bp_fdata, concated_level_data);
                            bpFileWriter.Close();
                            //std::cout << "processor " << rank << " finish FileWriter_ad\n";                            
                        }
                        free(concated_level_data);
                        count ++;
                        concated_level_size = 0;
                        prev_index = j + 1;
                    }
                }
                level_merged_count.push_back(merged_count);
            }
	    /*
	if(!rank){
        const std::map<std::string, adios2::Params> variables =
            bpIO.AvailableVariables();

        for (const auto variablePair : variables)
        {
            std::cout << "Name: " << variablePair.first;

            for (const auto &parameter : variablePair.second)
            {
                std::cout << "\t" << parameter.first << ": " << parameter.second
                          << "\n";
            }
        }
	}
	*/
            return level_merged_count;
        }

        void write_metadata(uint8_t const * metadata, uint32_t size) const {
            FILE * file = fopen(metadata_file.c_str(), "w");
            fwrite(metadata, 1, size, file);
            fclose(file);
        }

        ~HPSSFileWriter(){}

        void print() const {
            std::cout << "HPSS file writer." << std::endl;
        }
    private:

        uint32_t min_size = 0;
        std::vector<std::string> level_files;
        std::string metadata_file;
    };
}
#endif
