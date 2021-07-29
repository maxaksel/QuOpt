# QuOpt

### About
QuOpt is a parallelized C++ program designed to find optimal quantum circuits for user-defined unitary operations over three linearly connected qubits. Originally designed purely for the Toffoli gate, this software is useful for optimizing any frequently-occuring three-qubit operation.

Please direct any queries to Max Aksel Bowman (mbowman@anl.gov).

### Installation
QuOpt depends on MPI and Google Ceres (which itself relies on several software packages including but not limited to `Eigen`, `gflags`, and `glog`). `GTest` is also required for the test suite. Instructions on how to install Google Ceres can be found [here](http://ceres-solver.org/installation.html). Instructions on installing MPI can be found [here](https://www.mpich.org/downloads/). `GTest` must be built from [source](https://github.com/google/googletest). Many of these libraries utilize `cmake` to manage the build process, as does this project. Therefore, `cmake` must be installed first.

### Usage
The first step in using this program is editing program hyperparameters in `src/main.cpp`. These hyperparameters include number of parameters (and by extension number circuit of layers) and the number of Levenberg-Marquardt starting points to test.

`cd bin && mpirun -np 15 ./main` will run this program with MPI over 15 processes. Use this option for running over multiple cores on a single machine.

For running at a larger scale, please see the SLURM script in `bin`. This script is designed specifically to be run on the Bebop supercomputer at Argonne National Laboratory (`sbatch ./bebop_submit_slurm.sh`), but can be modified to suit the needs of any supercomputer using SLURM.

### Testing
To test this software, run the command `ctest` in the root directory.

### Findings
Coming soon!
