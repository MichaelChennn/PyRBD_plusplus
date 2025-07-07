## Recommended Environment

- WSL2 with Ubuntu 22.04
- Python 3.10 or higher
- GCC 11/G++ 11

## Install GCC, G++, CMake

```bash
sudo apt update
# sudo apt install -y gcc-11 g++-11 build-essential cmake
```

## Python Environment

```bash
conda create -n pyrbd_plusplus python=3.10
conda activate pyrbd_plusplus
pip install --upgrade pip
pip install -r requirements.txt
conda install -c conda-forge gcc_linux-64=11 gxx_linux-64=11

```

## Install Cpp Dependencies

```bash
rm -rf build
mkdir build
cd build
cmake -DCMAKE_C_COMPILER=x86_64-conda-linux-gnu-gcc -DCMAKE_CXX_COMPILER=x86_64-conda-linux-gnu-g++ ..
make distclean && make -j
```

## Insall Python Package

```bash
cd ..
pip install -e .
```
