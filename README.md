## Recommended Environment

- WSL2 with Ubuntu 22.04
- Python 3.10 or higher
- GCC 11/G++ 11

## Install GCC, G++, CMake

```bash
sudo apt update
```

## Python Environment

```bash
conda create -n pyrbd_plusplus python=3.10
conda activate pyrbd_plusplus
pip install --upgrade pip
pip install -r requirements.txt
conda install -c conda-forge gcc_linux-64=11 gxx_linux-64=11

```

## Install CPP Dependencies

```bash
chmod +x rebuild.sh
./rebuild.sh
```
