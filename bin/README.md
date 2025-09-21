# Software
Overview of codes used for systematic study of interchain and intrachain electrostatic repulsion. Other functions not used in the calculations can be found in the code; however, they have not been validated and may contain mistakes.

# Main files
Minimum working example in the above `demo` can be run by copying the contents of `./utils/` into the present working directory, running `$ make`, and then `$ ./pe-go`. 
The `Makefile` links the source C codes in `./src/`. The script `pe-go` is placed in the present working directory with the compiled codes in `./obj/`.

Polymer parameters backbone charge fraction $\alpha$ and polymer hydrophobicity $\chi$ can be set in `ac.txt`. 
`input.h` contains the remaining inputs needed for setting system size, bulk solution salt concentration, and convergence parameters. 
Physical constants, grid discretization lengths, and other unfrequently changed parameters are set in `./src/consts.h`.

While much of the scaffolding for multiblock polymers and multicomponent brushes have been set up in the code from previous projects, the calculations can only presently run for single-component, homogeneous PE brushes. 
Only small changes are needed in most of the code to consider these changes, but major changes are likely needed in the FFT section.

# Seed strategies
Although the `demo` does not require a previous run's results as an input, obtaining non-trivial PE brush morphologies often require seeds to ensure stable convergence. The parameter `init_opt` in `input.h` determines how the fields are initialized. 
- For `init_opt = 0`, the fields are seeded for a system with constant polymer fraction throughout the box. This is enough for very simple morphologies.
- `init_opt = 1` takes in a seed file of name specified in `input.h` (`Win[30]`) with the format of previously written `W_a*.dat` field files. Although these files written in 1D format, the dimensions of the original system are required to be inputted in `input.h`. If the current box is greater than the box of the seed, the edges of the previous calculations are used to fill in the areas of the larger box.
- `init_opt = -99` takes in a seed file but fills larger box sizes with fields corresponding to pure solvent.

# Convergence
For non-trivial PE brush morphologies, convergence using the techniques used in these codes can be very challenging. While Anderson Mixing is employed for accelerating convergence, several adaptations were made to increase stability. The mixing parameters `wand` and `wcmp` in `input.h` are the primary hyperparameters used to bias the algorithm towards speed or stability, respectively. `wand` generally determines how faithfully the calculations should follow the fields recommended by Anderson Mixing. `wcmp` modifies the fields based on how far each point is from satisfying the local incompressibility condition, $\phi_p + \phi_s = 1$. Generally, `wand = wcmp = 0.1` provides stable performance. `wand = 0.2, wcmp = 0.01` can be used when the seed file is thought to be relativelyclose to convergence.

# Batch runs
`./submit/` contains utilities for submitting single jobs (`submit.sh`, set up for use on Purdue Anvil CPU) and batch jobs (`submit_loop.c`). `submit_loop_in.txt` contains the changes between the batch jobs of the file on the first line (e.g., `input.h`). The integer preceeding `!` is the line number of the change and the text afterwards is what the line should be replaced with. The batch jobs are named after the parent folder and then the number of that copy (e.g., `folder.1`, `folder.2`). 

# Dependencies
The only setup necessary should be to ensure FFTW and OMP can be accessed by `make` through `lfftw3`, `fopenmp`, and `fftw3_omp`. For most computing clusters, this should be straightforward to implement given how common these packages are.
