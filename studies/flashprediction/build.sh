#!/bin/bash

# Build script for flash prediction executable

# Check if environment is set up
if [ -z "$LARCV_INCDIR" ]; then
    echo "Error: UBDL environment not set up!"
    echo "Please run:"
    echo "  cd /path/to/ubdl"
    echo "  source setenv_py3.sh"
    echo "  source configure.sh"
    exit 1
fi

# Create build directory if it doesn't exist
if [ ! -d "build" ]; then
    mkdir build
fi

cd build

# Configure with cmake
echo "Configuring with cmake..."
cmake ..

# Build
echo "Building..."
make -j$(nproc)

# Check if build was successful
if [ $? -eq 0 ]; then
    echo "Build successful!"
    echo "Executables created:"
    echo "  - build/flashprediction/calculate_flash_predictions"
    echo "  - build/flashprediction/run_flashprediction"
    echo ""
    echo "Test the executable:"
    echo "  ./build/flashprediction/calculate_flash_predictions -h"
else
    echo "Build failed!"
    exit 1
fi