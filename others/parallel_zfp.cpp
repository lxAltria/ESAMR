#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <zfp.h>
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

unsigned char * zfp_compress_3D(float * array, double tolerance, size_t r1, size_t r2, size_t r3, size_t *out_size){
	int status = 0;    /* return value: 0 = success */
	zfp_type type;     /* array scalar type */
	zfp_field* field;  /* array meta data */
	zfp_stream* zfp;   /* compressed stream */
	void* buffer;      /* storage for compressed stream */
	size_t bufsize;    /* byte size of compressed buffer */
	bitstream* stream; /* bit stream to write to or read from */
	size_t zfpsize;    /* byte size of compressed stream */

	/* allocate meta data for the 3D array a[nz][ny][nx] */
	type = zfp_type_float;
	field = zfp_field_3d(array, type, r3, r2, r1);

	/* allocate meta data for a compressed stream */
	zfp = zfp_stream_open(NULL);

	/* set compression mode and parameters via one of three functions */
	/*  zfp_stream_set_rate(zfp, rate, type, 3, 0); */
	/*  zfp_stream_set_precision(zfp, precision); */
	zfp_stream_set_accuracy(zfp, tolerance);

	/* allocate buffer for compressed data */
	bufsize = zfp_stream_maximum_size(zfp, field);
	buffer = malloc(bufsize);

	/* associate bit stream with allocated buffer */
	stream = stream_open(buffer, bufsize);
	zfp_stream_set_bit_stream(zfp, stream);
	zfp_stream_rewind(zfp);

	zfpsize = zfp_compress(zfp, field);

	zfp_field_free(field);
	zfp_stream_close(zfp);
	stream_close(stream);
	*out_size = zfpsize;
	return (unsigned char *)buffer;
}

float * zfp_decompress_3D(unsigned char * comp_data, double tolerance, size_t buffer_size, size_t r1, size_t r2, size_t r3){
	int status = 0;    /* return value: 0 = success */
	zfp_type type;     /* array scalar type */
	zfp_field* field;  /* array meta data */
	zfp_stream* zfp;   /* compressed stream */
	void* buffer;      /* storage for compressed stream */
	size_t bufsize;    /* byte size of compressed buffer */
	bitstream* stream; /* bit stream to write to or read from */
	size_t zfpsize;    /* byte size of compressed stream */

	/* allocate meta data for the 3D array a[nz][ny][nx] */
	float * array = (float *) malloc(r1 * r2 * r3 * sizeof(float));
	type = zfp_type_float;
	field = zfp_field_3d(array, type, r3, r2, r1);

	/* allocate meta data for a compressed stream */
	zfp = zfp_stream_open(NULL);

	/* set compression mode and parameters via one of three functions */
	/*  zfp_stream_set_rate(zfp, rate, type, 3, 0); */
	/*  zfp_stream_set_precision(zfp, precision); */
	zfp_stream_set_accuracy(zfp, tolerance);

	/* allocate buffer for compressed data */
	bufsize = zfp_stream_maximum_size(zfp, field);
	// buffer = malloc(bufsize);
	buffer = (void *) comp_data;
	bufsize = buffer_size;

	/* associate bit stream with allocated buffer */
	stream = stream_open(buffer, bufsize);
	zfp_stream_set_bit_stream(zfp, stream);
	zfp_stream_rewind(zfp);

	zfp_decompress(zfp, field);

	zfp_field_free(field);
	zfp_stream_close(zfp);
	stream_close(stream);
	return array;
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
        adios2::IO bpIO = ad.DeclareIO("WriteZFP");
        std::string zip_filename = "p" + std::to_string(size) + "_zfp/eb_1e-" + std::to_string(i+1);
	    // Compress Input Data
	    if (rank == 0){
            printf("Compressing eb: %f\n", abs_bound[i]);
        }
	    MPI_Barrier(MPI_COMM_WORLD);
	    if(rank == 0){
            compress_time[i] = - MPI_Wtime();
        }
        // true compression API
	    unsigned char * compressed_data = zfp_compress_3D(data, abs_bound[i], n1, n2, n3, &compressed_size[i]); 
		MPI_Barrier(MPI_COMM_WORLD);
    	if(rank == 0){
	    	compress_time[i] += MPI_Wtime();
	    }
	    MPI_Barrier(MPI_COMM_WORLD);
	    if(rank == 0){
            write_zip_time[i] = - MPI_Wtime();
        }
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
        free(compressed_data);
    }
    free(data);
	for(int i=0; i<num_reb; i++){
        adios2::ADIOS ad(MPI_COMM_WORLD);
        adios2::IO readIO = ad.DeclareIO("ReadZFP");
        std::string zip_filename = "p" + std::to_string(size) + "_zfp/eb_1e-" + std::to_string(i+1);
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
	    float * decompressed_data = zfp_decompress_3D(compressed_data, abs_bound[i], compressed_size[i], n1, n2, n3); 
	    MPI_Barrier(MPI_COMM_WORLD);
	    if(rank == 0){
		    decompress_time[i] += MPI_Wtime();
	    }
        free(compressed_data);
        free(decompressed_data);
    }

	if (rank == 0){
		printf ("ZFP Finish parallel compressing\n");
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
	MPI_Finalize();
	return 0;
}
