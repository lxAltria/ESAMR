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
mkdir -p refactored_data<br />
Refactor: ./test/test_refactor $data_file $data_type $encode_opt $reorder_opt $num_dims $dim0 $dim1 $dim2<br />
./test/test_refactor ../external/SZ3/data/Uf48.bin.dat 0 3 0 1 3 100 500 500<br />
Retrieval: ./test/test_retrieval $data_file $data_type $error_mode $error $dims $dim0 $dim1 $dim2<br />
./test/test_retrieval ../external/SZ3/data/Uf48.bin.dat 0 1 1.0 3 100 500 500<br />