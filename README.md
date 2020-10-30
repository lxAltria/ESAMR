# Multi-precision data refactoring framework
Project: ESAMR (Enabling Scalable Analysis using Multi-precision Refactoring)<br />
Sponsor: ORNL LDRD<br />
Authors: Xin Liang, Jieyang Chen, Lipeng Wan, Ben Whitney, Qian Gong<br />
Supervisors: Scott Klasky, Rick Archibald<br />

# Installation
git clone https://github.com/lxAltria/Multiprecision-data-refactoring.git<br />
cd Multiprecision-data-refactoring<br />
./build_script.sh<br />

# Usage
cd build<br />
mkdir -p refactor_data<br />
Refactor: ./test/test_refactor $data_file $data_type $encode_opt $reorder_opt $num_dims $dim0 $dim1 $dim2<br />
./test/test_refactor ../external/SZ3/data/Uf48.bin.dat 0 3 0 1 3 100 500 500<br />
Retrieval: ./test/test_retrieval $data_file $data_type $error_mode $error $dims $dim0 $dim1 $dim2<br />
./test/test_retrieval ../external/SZ3/data/Uf48.bin.dat 0 1 1.0 3 100 500 500<br />

# Notes and Parameters
During refactoring, the location of refactored data is hardcoded to "refactor_data/" directory under current directory. Need to create the directory before writing.<br />
During retrieving, the location of recomposed data is hardcded to "mgard.recompose" under current directory.<br />
data_file: path to input date file.<br />
data_type: 0 for float, 1 for double.<br />
encode_opt: encode options for bitplanes (see include/data_enc.hpp)<br />
0: default bitplane encoding<br />
1: default bitplane encoding with sign postpone, i.e. encoding the sign after first 1 of current data<br />
2: runlength encoding<br />
3: mix of default bitplane encoding and runlength encoding<br />
4: mix of default bitplane encoding with sign postpone and runlength encoding<br />
reorder_opt: reordering options for bitplanes (see include/error_est.hpp)<br />
0: in order, which records bitplanes following level ordering (the same as MGARD)<br />
1: round-robin, which iteratively picks one bitplane for each level<br />
2: uniform error, which picks bitplane such that every level will have the same L-infinity error<br />
3: greedy shuffling, which picks bitplane according to a greedy method<br />
error mode: error metric during retreival (see include/error_est.hpp)<br />
1: max error, i.e. L-infty<br />
2: squared error, i.e. L-2<br />
3: PSNR<br />
