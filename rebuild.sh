#!/bin/bash
# quick_rebuild.sh - Quick rebuild script

echo "🚀 Quick rebuild for PyRBD++..."

# Check if we're in the right directory
if [ ! -f "setup.py" ]; then
    echo "❌ Please run from project root directory (containing setup.py)"
    exit 1
fi

# Save starting directory
START_DIR=$(pwd)

# Execute build steps
echo "📁 Current directory: $START_DIR"
rm -rf build/ && \
mkdir build && \
cd build && \
cmake -DCMAKE_C_COMPILER=x86_64-conda-linux-gnu-gcc -DCMAKE_CXX_COMPILER=x86_64-conda-linux-gnu-g++ .. && \
make distclean && \
make -j$(nproc) && \
cd "$START_DIR" && \
pip install -e . && \
python -c "import pyrbd_plusplus; print('✅ Build successful!')" || echo "❌ Build failed"

echo "📍 Final directory: $(pwd)"