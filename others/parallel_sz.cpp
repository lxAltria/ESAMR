#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sz.h>
#include <mpi.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstring>
#include <vector>
#include <adios2.h>

template<class T>
void posix_read(std::string filename, T* data, size_t num_elements){
    int fd = open(filename.c_str(), O_RDONLY);
    read(fd, data, num_elements*sizeof(T));
    close(fd);
}

int main(int argc, char ** argv){
	MPI_Init(&argc, &argv);
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);	
		
    if (rank == 0) printf ("Start parallel compressing using %d cores... \n", size);

    // dataset information
    std::string filename = argv[1];
    const int n1 = atoi(argv[2]);
    const int n2 = atoi(argv[3]);
    const int n3 = atoi(argv[4]);
    int num_elements = n1 * n2 * n3;
    const int num_reb = 7;
    double rel_bound[num_reb] = {0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001};

    // timing
	double start, end;
    double read_ori_file_time = 0;
    double write_zip_time[num_reb] = {0};
    double read_zip_time[num_reb] = {0};
    double compress_time[num_reb] = {0};
    double decompress_time[num_reb] = {0};

	size_t compressed_size[num_reb] = {0};

	int status = 0;
	float * data = NULL;
	SZ_Init(NULL);
    // Read Input Data
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0){
        read_ori_file_time = - MPI_Wtime();
        data = (float*)malloc(num_elements* sizeof(float));             
        posix_read(filename, data, num_elements);           
        MPI_Bcast(&num_elements, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(data, num_elements, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }
    else{
        MPI_Bcast(&num_elements, 1, MPI_INT, 0, MPI_COMM_WORLD);
        data = (float *) malloc(num_elements * sizeof(float));
        MPI_Bcast(data, num_elements, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0){
        read_ori_file_time += MPI_Wtime();
    }
    // compute abs_bound
    float max_v = data[0];
    float min_v = data[0];
    for(int i=1; i<num_elements; i++){
	    if(data[i] > max_v) max_v = data[i];
	    if(data[i] < min_v) min_v = data[i];
    }
    float value_range = max_v - min_v;
    std::vector<double> abs_bound(num_reb);
    for(int i=0; i<num_reb; i++){
	    abs_bound[i] = value_range * rel_bound[i];
    }
    for (int i=0; i<num_reb; i++) {
        adios2::ADIOS ad(MPI_COMM_WORLD);
        adios2::IO bpIO = ad.DeclareIO("WriteSZ");
        std::string zip_filename = "p" + std::to_string(size) + "_sz/eb_1e-" + std::to_string(i+1);
	    // Compress Input Data
	    if (rank == 0){
            printf("Compressing using absolute eb: %f\n", abs_bound[i]);
        }
	    MPI_Barrier(MPI_COMM_WORLD);
	    if(rank == 0){
            compress_time[i] = - MPI_Wtime();
        }
        // true compression API
	    unsigned char * compressed_data = SZ_compress_args(SZ_FLOAT, data, &compressed_size[i], ABS, abs_bound[i], abs_bound[i], abs_bound[i], 0, 0, n1, n2, n3); 
		MPI_Barrier(MPI_COMM_WORLD);
    	if(rank == 0){
	    	compress_time[i] += MPI_Wtime();
	    }
	    MPI_Barrier(MPI_COMM_WORLD);
	    if(rank == 0){
            write_zip_time[i] = - MPI_Wtime();
        }
        // adios write 
        {
            uint64_t shape = compressed_size[i] * size;
            uint64_t start = compressed_size[i] * rank;
            uint64_t length = compressed_size[i];
            adios2::Variable<uint8_t> bp_fdata = bpIO.DefineVariable<uint8_t>(zip_filename, {shape}, {start}, {length}, adios2::ConstantDims);
            adios2::Engine bpFileWriter = bpIO.Open(zip_filename, adios2::Mode::Write);
            bpFileWriter.Put<uint8_t>(bp_fdata, compressed_data);
            bpFileWriter.Close();
        }
	    MPI_Barrier(MPI_COMM_WORLD);
	    if(rank == 0){
		    write_zip_time[i] += MPI_Wtime();
	    }
	    if(rank == 0){
		    printf("writing to %s\n", zip_filename.c_str());
	    }
        free(compressed_data);
    }
    free(data);
	for(int i=0; i<num_reb; i++){
        adios2::ADIOS ad(MPI_COMM_WORLD);
        adios2::IO readIO = ad.DeclareIO("ReadSZ");
        std::string zip_filename = "p" + std::to_string(size) + "_sz/eb_1e-" + std::to_string(i+1);
        unsigned char * compressed_data = (unsigned char *) malloc(compressed_size[i]);
        // Read compressed Data
        MPI_Barrier(MPI_COMM_WORLD);
        if(rank == 0){
            read_zip_time[i] = - MPI_Wtime();
        }
        {
            adios2::Engine bpFileReader = readIO.Open(zip_filename, adios2::Mode::Read);
            adios2::Variable<uint8_t> bp_fdata = readIO.InquireVariable<uint8_t>(zip_filename);
            uint64_t start = compressed_size[i] * rank;
            uint64_t length = compressed_size[i];
            bp_fdata.SetSelection(adios2::Box<adios2::Dims>({start}, {length}));
            bpFileReader.Get<uint8_t>(bp_fdata, compressed_data, adios2::Mode::Sync);
            bpFileReader.Close();
        }        
        MPI_Barrier(MPI_COMM_WORLD);
        if(rank == 0){
            read_zip_time[i] += MPI_Wtime();
        }
    	// Decompress Compressed Data
	    MPI_Barrier(MPI_COMM_WORLD);
	    if(rank == 0){
            decompress_time[i] = - MPI_Wtime();
        }
	    float * decompressed_data = (float*)SZ_decompress(SZ_FLOAT, compressed_data, compressed_size[i], 0, 0, n1, n2, n3); 
	    MPI_Barrier(MPI_COMM_WORLD);
	    if(rank == 0){
		    decompress_time[i] += MPI_Wtime();
	    }
        free(compressed_data);
        free(decompressed_data);
    }

	if (rank == 0){
		printf ("SZ Finish parallel compressing\n");
        printf ("Timecost of reading input file = %.4f\n", read_ori_file_time);
        uint32_t accum_size = 0;
        for (int i=0; i<num_reb; i++) {
            accum_size += compressed_size[i];
            printf("eb = %.7f\n", rel_bound[i]);
            printf("Compression ratio = %.4f, accumulated ratio = %.4f\n", 1.0*num_elements*sizeof(float) / compressed_size[i], 1.0*num_elements*sizeof(float) / accum_size);
            printf("Timecost of compressing using %d processes = %.6f seconds\n", size, compress_time[i]);
            printf("Timecost of decompressing using %d processes = %.6f seconds\n", size, decompress_time[i]);
            printf("Timecost of reading compressed files = %.6f seconds\n", read_zip_time[i]);
            printf("Timecost of writing compressed files = %.6f seconds\n\n", write_zip_time[i]);
		}
	}
	SZ_Finalize();
	MPI_Finalize();
	return 0;
}
