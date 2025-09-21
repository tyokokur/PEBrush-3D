# Software
Overview of codes used for systematic study of interchain and intrachain electrostatic repulsion. Other functions not used in the calculations can be found in the code; however, they have not been validated and may contain mistakes.

# Main files
Minimum working example in `demo` can be run by copying the contents of `./utils/` into the present working directory, running `$ make`, and then `$ ./pe-go`. 
Polymer parameters backbone charge fraction $\alpha$ and polymer hydrophobicity $\chi$ can be changed in `ac.txt`. 
`input.h` contains the remaining inputs needed for setting system size, bulk solution salt concentration, and convergence parameters. 
Physical constants, grid discretization lengths, and other unfrequently changed parameters are set in `./src/consts.h`.

# Setup
The only setup necessary should be to specify the location of the FFTW library on your setup. 
