#!/bin/sh

current_dir=`pwd`
hts_ld= echo $LD_LIBRARY_PATH | grep -o $current_dir/htslib/
clang_prefix=/install

clang_ld= echo $LD_LIBRARY_PATH | grep -o $clang_prefix/lib
clang_path= echo $PATH | grep -o $clang_prefix/bin

if [ "$hts_ld" != "$current_dir/hstlib/" ]
then
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$current_dir/htslib/
fi

if [ "$clang_ld" != "$clang_prefix/lib" ]
then
    export LD_LIBRARY_PATH=$clang_prefix/lib:$LD_LIBRARY_PATH
fi

if [ "$clang_path" != "$clang_prefix/bin" ]
then
    export PATH=$clang_prefix/bin:$PATH
fi
