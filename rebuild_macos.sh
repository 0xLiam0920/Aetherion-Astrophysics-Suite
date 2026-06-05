#!/bin/bash

# Navigate to the project root directory
cd "$HOME/aetherionsuite" || exit 1

# Ensure submodules are populated
git submodule update --init --recursive

# Generate build files and compile using 4 CPU cores
cmake -S . -B build && cmake --build build -j4

# Navigate to the build output directory
cd build || exit 1

# Run the application (checks for Mac bundle first, falls back to binary)
if [ -d "blackhole-sim.app" ]; then
    open blackhole-sim.app
else
    ./blackhole-sim
fi

