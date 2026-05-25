## I hated typing three commands (yes im lazy) so I made this script. Run it wherever you want.

# cd's in ~/aetherionsuite and runs the build script, then the build directory and ./blackhole-sim
cd ~/aetherionsuite
cmake -S . -B build && make -C build -j4 ## if you have a different number of CPU cores, change the 4 to something else for faster builds. 4 is fine for now
cd build
./blackhole-sim
# note: if for some weird reason this script isn't executable, try running `chmod +x rebuld_macos.sh` in the terminal. You only need to do that once.