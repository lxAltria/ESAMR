#ifndef _MDR_HPSS_FILE_RETRIEVER_HPP
#define _MDR_HPSS_FILE_RETRIEVER_HPP

#include "RetrieverInterface.hpp"
#include <cstdio>
#include "RefactorUtils.hpp"
#include <mpi.h>
#include <adios2.h>

namespace MDR {
    // Data retriever for files
    class HPSSFileRetriever : public concepts::RetrieverInterface {
    public:
        HPSSFileRetriever(const std::string& metadata_file, const std::vector<std::string>& level_files) : metadata_file(metadata_file), level_files(level_files) {
            bitplane_offsets = std::vector<uint32_t>(level_files.size(), 0);
            segment_offsets = std::vector<uint32_t>(level_files.size(), 0);
        }

        void init(const std::vector<std::vector<uint32_t>>& level_sizes, const std::vector<std::vector<double>>& level_squared_errors, const std::vector<std::vector<uint32_t>>& level_merged_count){
            segment_bitplane_count = level_merged_count;
            level_segment_size.clear();
            level_segment_error.clear();
            for(int i=0; i<level_sizes.size(); i++){
                std::vector<uint32_t> segment_size;
                std::vector<double> segment_error;
                segment_error.push_back(level_squared_errors[i][0]);
                uint32_t bitplane_index = 0;
                for(int j=0; j<level_merged_count[i].size(); j++){
                    uint64_t size = 0;                    
                    for(int k=0; k<level_merged_count[i][j]; k++){
                        size += level_sizes[i][bitplane_index + k];
                    }
                    bitplane_index += level_merged_count[i][j];
                    if(size > (1ULL << 32)){
                        std::cerr << "Size exceed 32 bit range\n";
                        exit(-1);
                    }
                    segment_size.push_back(size);
                    segment_error.push_back(level_squared_errors[i][bitplane_index]);
                }
                level_segment_size.push_back(segment_size);
                level_segment_error.push_back(segment_error);
            }
        }

        std::vector<std::vector<const uint8_t*>> retrieve_level_components(const std::vector<std::vector<uint32_t>>& level_sizes, const std::vector<uint32_t>& retrieve_sizes, const std::vector<uint8_t>& prev_level_num_bitplanes, const std::vector<uint8_t>& level_num_bitplanes){
            int rank, size;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_Comm_size(MPI_COMM_WORLD, &size);
            adios2::ADIOS ad(MPI_COMM_WORLD);
            adios2::IO readIO = ad.DeclareIO("ReadPrecisionSegments");

            std::vector<std::vector<const uint8_t*>> level_components;
            release();
            uint32_t retrieve_size = 0;
            for(int i=0; i<level_files.size(); i++){
                // std::cout << "Retrieve " << +level_num_bitplanes[i] << " (" << +(level_num_bitplanes[i] - prev_level_num_bitplanes[i]) << " more) bitplanes from level " << i << std::endl;
                std::vector<const uint8_t*> interleaved_level;
                for(int j=prev_level_num_bitplanes[i]; j<level_num_bitplanes[i]; j++){
                    // read one segment 
                    uint32_t num_bitplanes = segment_bitplane_count[i][j];
                    uint8_t * buffer = (uint8_t *) malloc(level_segment_size[i][segment_offsets[i]]);
                    // FILE * file = fopen((level_files[i] + "_" + std::to_string(segment_offsets[i])).c_str(), "r");
                    // fread(buffer, sizeof(uint8_t), level_segment_size[i][segment_offsets[i]], file);
                    // fclose(file);
                    std::string filename = level_files[i] + "_" + std::to_string(segment_offsets[i]);
                    // MPIIO_read(MPI_COMM_WORLD, rank, size, level_segment_size[i][segment_offsets[i]], buffer, filename);
                    {
            			//adios2::IO readIO = ad.DeclareIO(filename);
            			MPI_Barrier(MPI_COMM_WORLD);
                        adios2::Engine bpFileReader = readIO.Open(filename, adios2::Mode::Read);
			/*
			if(bpFileReader){
				if(!rank){
					const std::map<std::string, adios2::Params> variables =
            				readIO.AvailableVariables();
        				for (const auto variablePair : variables)
        				{
            					std::cout << "Name: " << variablePair.first;

            					for (const auto &parameter : variablePair.second)
            					{
                					std::cout << "\t" << parameter.first << ": " << parameter.second << "\n";
            					}
        				}

				}
			}
			else{
				std::cout << "open file failed\n";
				std::cout << rank << " " << size << std::endl;
				MPI_Abort(MPI_COMM_WORLD, -1);
			}
			*/
                        adios2::Variable<uint8_t> bp_fdata = readIO.InquireVariable<uint8_t>(filename);
			
			if(bp_fdata){
				std::cout << filename << "_found in rank " << rank << "\n";
				fflush(stdout);
			}
			else{
				std::cout << filename << "_not_found in rank " << rank << std::endl;
				fflush(stdout);
					const std::map<std::string, adios2::Params> variables =
                                        readIO.AvailableVariables();
                                        for (const auto variablePair : variables)
                                        {
                                                std::cout << "Name: " << variablePair.first;

                                                for (const auto &parameter : variablePair.second)
                                                {
                                                        std::cout << "\t" << parameter.first << ": " << parameter.second << "\n";
                                                }
                                        }
				std::cout << "list of existing variables done\n";
				MPI_Abort(MPI_COMM_WORLD, 0);
			}
			
                        bp_fdata.SetSelection(adios2::Box<adios2::Dims>({(uint64_t)level_segment_size[i][segment_offsets[i]] * rank}, {(uint64_t)level_segment_size[i][segment_offsets[i]]}));
                        bpFileReader.Get<uint8_t>(bp_fdata, buffer, adios2::Mode::Sync);
                        bpFileReader.Close();
                    }
                    concated_level_components.push_back(buffer);
                    // interleave level component
                    const uint8_t * pos = buffer;
                    for(int k=0; k<num_bitplanes; k++){
                        interleaved_level.push_back(pos);
                        pos += level_sizes[i][bitplane_offsets[i] + k];
                    }
                    bitplane_offsets[i] += num_bitplanes;
                    retrieve_size += level_segment_size[i][segment_offsets[i]];
                    segment_offsets[i] ++;
                }
                level_components.push_back(interleaved_level);
            }
            uint32_t total_retrieve_size = 0;
            for(int i=0; i<level_files.size(); i++){
                for(int j=0; j<level_num_bitplanes[i]; j++){
                    total_retrieve_size += level_segment_size[i][j];
                }
            }
            // std::cout << "Total retrieve size = " << total_retrieve_size << std::endl;
            return level_components;
        }

        uint8_t * load_metadata() const {
            FILE * file = fopen(metadata_file.c_str(), "r");
            fseek(file, 0, SEEK_END);
            uint32_t num_bytes = ftell(file);
            rewind(file);
            uint8_t * metadata = (uint8_t *) malloc(num_bytes);
            fread(metadata, 1, num_bytes, file);
            fclose(file);
            return metadata;
        }

        void release(){
            for(int i=0; i<concated_level_components.size(); i++){
                free(concated_level_components[i]);
            }
            concated_level_components.clear();
        }

        const std::vector<std::vector<uint32_t>>& get_level_segment_size(){
            return level_segment_size;
        }

        const std::vector<std::vector<double>>& get_level_segment_error(){
            return level_segment_error;
        }

        ~HPSSFileRetriever(){}

        void print() const {
            std::cout << "Merged file retriever." << std::endl;
        }
    private:
        std::vector<std::string> level_files;
        std::string metadata_file;
        std::vector<uint32_t> bitplane_offsets; // tracking the next bitplane index for each level
        std::vector<uint32_t> segment_offsets; // tracking the next segment index for each level
        std::vector<uint8_t*> concated_level_components;

        std::vector<std::vector<uint32_t>> segment_bitplane_count;
        std::vector<std::vector<uint32_t>> level_segment_size;
        std::vector<std::vector<double>> level_segment_error;
    };
}
#endif
