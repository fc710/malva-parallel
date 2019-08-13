# MALVA

Parallel version of Malva tool found at https://github.com/AlgoLab/malva

## Prerequisite

###  NVIDIA GPU

A CUDA-enabled NVIDIA GPU is required to be able to build and run the application.
If you are not sure which GPU is installed in your system, run :
```bash
lspci  -v -s  $(lspci | grep ' VGA ' | cut -d" " -f 1)
```
You also need to know its compute capability: a list with the information can be found [here](https://developer.nvidia.com/cuda-gpus).
For example the GeForce GTX 1080 has a compute capability 6.1.

## Building the Compiler

### 0. Compiler with OpenMP offloading to NVIDIA GPUs

Compilers pre-installed on the system(such as gcc) or downloaded pre-built binaries support OpenMP by default, but don't support OpenMP offloading to other devices, so a custom version needs to be built.
Currently **gcc** doesn't support [OpenMP nvptx offloading](https://gcc.gnu.org/wiki/Offloading#nvptx_offloading), so we'll use **clang** from the [LLVM project](https://llvm.org/).

### 1. Build Prerequisites

Before building LLVM make sure the following software is installed:

- Standard tools *make*, *cmake*, *tar*, *xz*.
- An already working compiler such as *gcc*. Regarding the version, the latest *gcc 9*, can cause issues when building the OpenMP runtime, and the latest **working** version of the CUDA Toolkit (10.0) supports up to *gcc 7*, so an installation of *gcc 7* is recommended for the building process.
- For the runtime libraries *libelf*.
- The NVIDIA CUDA toolkit. It is recommended the [version 10.0](https://developer.nvidia.com/cuda-10.0-download-archive) , since *clang 8.0.1* doesn't seem to support version 10.1.

### 2. Get the sources

Download the LLVM sources for the current latest version [8.0.1](https://releases.llvm.org/download.html#8.0.1) or for other versions
(not recommended) from https://releases.llvm.org/.
Several source components are available, the **required** are:
- LLVM source code
- Clang source code
- compiler-rt source code
- OpenMP source code

Extract the tarballs. This should result in 4 directories: *llvm-8.0.1.src*, *cfe-8.0.1.src*, *compiler-rt-8.0.1.src*, *openmp-8.0.1.src*. Run the following commands to allow in-tree building:

```bash
mv cfe-8.0.1.src llvm-8.0.1.src/tools/clang
mv compiler-rt-8.0.1.src llvm-8.0.1.src/projects/compiler-rt
mv openmp-8.0.1.src llvm-8.0.1.src/projects/openmp
```

### 3. Build the Compiler

Choose the location where the compiler will be installed (in this example */install*).

```bash
install_prefix=/install
```

Recall the compute capability of your GPU (e.g for *6.1* you'll use sm_61), and of all GPUs you want to support.

Now you're ready to build:

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$install_prefix -DCLANG_OPENMP_NVPTX_DEFAULT_ARCH=sm_61 -DLIBOMPTARGET_NVPTX_COMPUTE_CAPABILITIES=35,60,61,70 ../llvm-8.0.1.src
```
#### Note:
*CMAKE_INSTALL_PREFIX* specifies where the final binaries and libraries will be installed.
*CLANG_OPENMP_NVPTX_DEFAULT_ARCH* should use your GPU compute capability. *LIBOMPTARGET_NVPTX_COMPUTE_CAPABILITIES* all the ones you want your compiler to support.

Now the slow part: if your CPU has multiple cores, you can use the flag -j<*number of cores*> to speed it up.

```bash
make -j8
make -j8 install
```

### 4. Rebuild the OpenMP Runtime

You'll need to rebuild the OpenMP project, using Clang previously built - or alternatively repeat all the building process of step 3 (bootstrapping Clang). To just rebuild OpenMP:

```bash
cd ..
mkdir build-openmp
cd build-openmp
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$install_prefix-DCMAKE_C_COMPILER=/install/bin/clang -DCMAKE_CXX_COMPILER=/install/bin/clang++ -DCMAKE_C_FLAGS='--cuda-path=<*path/to/cuda*> --gcc-toolchain=<*/gcc/prefix*>' -DCMAKE_CXX_FLAGS='--cuda-path=/opt/cuda-10.0 --gcc-toolchain=<*/gcc/prefix*>' -DCUDA_TOOLKIT_ROOT_DIR=<*/path/to/cuda*> -DCUDA_HOST_COMPILER=<*path/to/gcc*> -DLIBOMPTARGET_NVPTX_COMPUTE_CAPABILITIES=35,60,61,70 ../llvm-8.0.1.src/projects/openmp
```

### Note carefully:
*CMAKE_C_FLAGS*, *CMAKE_CXX_FLAGS*, *CUDA_TOOLKIT_ROOT_DIR*, *CUDA_HOST_COMPILER* are optional, but necessary if there are version errors between clang and gcc/cuda.

- *CMAKE_C/CXX_COMPILER* contains the path to the newly built clang binaries.
- *CMAKE_C/CXX_FLAGS* is used to fix breakage if unsupported gcc libraries(version 9+) are selected by clang and to specify where to locate CUDA if it cannot be detected automatically:
  - *--cuda-path* specifies the path the CUDA toolkit.
  - *--gcc-toolchain* requires the prefix to the gcc toolchain, usually /usr if the gcc binary is in /usr/bin/. The version cannot be specified, so if multiple versions are installed in the same prefix, a fake install is needed.
  For example if i want to select gcc 8.3 libraries, assuming they are in /usr :
  
  ```bash
  version=8.3.0
  sudo mkdir -p $install_prefix/gcc/$version/include/c++
  sudo mkdir -p $install_prefix/gcc/$version/lib/gcc/x86_64-unknown-linux-gnu
  sudo ln -s /usr/include/c++/$version $install_prefix/gcc/$version/include/c++/$version
  sudo ln -s /usr/lib/gcc/x86_64-linux-gnu/$version $install_prefix/gcc/$version/lib/gcc/x86_64-unknown-linux-gnu/$version
  ```
  this way the flag is --gcc-toolchain=$install_prefix/gcc/$version
  - *CUDA_HOST_COMPILER* specifies the binary to a gcc 7 or lower for CUDA 10.0.

#### Build and install the OpenMP runtime libraries

```bash
make -j8
make -j8 install
```
## Download and Install Malva

#### Download the source code 

```bash
git clone --single-branch --branch gpu-accel --recursive https://github.com/fc710/malva-parallel.git
cd malva-parallel
```

#### Install the 3rd party libraries

```bash
cd sdsl-lite/build
./build.sh
cd ../../KMC
make
cd ../htslib
make
cd ..
```

#### Setup and run the config script to fix the environment variables

First tell the script where you installed clang by executing:
```bash
sed -i "5iclang_prefix=$install_prefix" config.sh
```
Then run :
```bash
./config.sh
```
Unless these environment variables are set permanentely (e.g added in ~/.bashrc), this script has to be rerun every time the variables are resetted.

#### Build the application

```bash
make
```

## Usage

Some examples on the usage [here](https://github.com/AlgoLab/malva#usage).