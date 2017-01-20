spatiocyte-cuda-zorder
======================

A CUDA standalone Spatiocyte package with Z-order curve voxels.

1. Get spatiocyte-cuda-zorder
    * $ git clone https://github.com/satya-arjunan/spatiocyte-cuda-zorder.git

2. Compile spatiocyte-cuda
    * add CUDA nvcc path (e.g., /usr/local/cuda-8.0/bin) into the PATH enviroment variable
    * $ make (or 'make -jx' with x the number of available CPU cores)

3. Run spatiocyte-cuda
    * $ ./spatiocyte-cuda

4. Visualize logged data
    * $ ./visualizer
