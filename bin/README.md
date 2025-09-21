# Software
Overview of codes used for systematic study of interchain and intrachain electrostatic repulsion. Other functions not used in the calculations can be found in the code; however, they have not been validated and may contain mistakes.

# Main files
Minimum working example in the above `demo` can be run by copying the contents of `./utils/` into the present working directory, running `$ make`, and then `$ ./pe-go`. 
Polymer parameters backbone charge fraction $\alpha$ and polymer hydrophobicity $\chi$ can be changed in `ac.txt`. 
`input.h` contains the remaining inputs needed for setting system size, bulk solution salt concentration, and convergence parameters. 
Physical constants, grid discretization lengths, and other unfrequently changed parameters are set in `./src/consts.h`.

# Batch runs
`./submit/` contains utilities for submitting single jobs (`submit.sh`, set up for use on Purdue Anvil CPU) and batch jobs (`submit_loop.c`). `submit_loop_in.txt` contains the changes between the batch jobs of the file on the first line (e.g., `input.h`). The integer preceeding `!` is the line number of the change and the text afterwards is what the line should be replaced with. The batch jobs are named after the parent folder and then the number of that copy (e.g., `folder.1`, `folder.2`). 

# Dependencies
The only setup necessary should be to ensure FFTW and OMP can be accessed by `make` through `lfftw3`, `fopenmp`, and `fftw3_omp`. For most computing clusters, this should be straightforward to implement given how common these packages are.
