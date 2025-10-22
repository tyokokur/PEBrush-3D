# PEBrush-3D
Codes and associated visualization codes for full (3D), self-consistent field theory (SCFT) calculations for polyelectrolyte brushes. C codes provide equilibrium polymer density distributions, solvent density distributions, ionic distributions, and electrostatic potential. Can be straightforwardly adapted for the multiblock charged macromolecular model used to described proteins. Corresponding equations can be found in `theory.pdf`.

Author: Yokokura, T. J. as part of the Prof. Rui Wang Group in the Chemical and Biomolecular Engineering Dept. at UC Berkeley. Please direct questions at ruiwang325(at)berkeley.edu.

# Sample Results 
Once converged, the calculations yield polymer density distributions for the system at equilibrium. For example, the effect of grafting density can be explored. Sparse brushes with hydrophobic interactions dominating tend to form clusters of pinned micelles (left); dense brushes stretch far into the solution (right).

<img src=/vis/sample-micelle.png height="225"/> <img src=/vis/sample-homogeneous.png height="300"/> 

Lateral inhomogeneities in only one axis of the grafting plane can also be explored. For example, cylinder-like aggregates formed by highly charged polyelectrolytes (left) can further melt to form stripes with tops deformed into beads (right).

<img src=/vis/sample-cylinder.png height="250"/> <img src=/vis/sample-cylstripe.png height="275"/> 

# File tree
- `bin/`: C codes used for SCFT calculations
- `demo/`: Example results generated for a systematic study of PE brushes under varying interchain and intrachain electrostatic repulsion
- `LICENCE`: Standard MIT license
- `readme.md`: This readme file
- `theory-full.pdf`: Derivation of self-consistent field equations
- `theory.pdf`: Concise derivation of key equations following standard buildup (see G. H. Fredrickson, Oxford Press, 2006)
- `vis/`: Python codes used for visualizations
  
# Overview
_Disclaimer: I am a chemical engineer, not a coder_

Dependencies
- SCFT calculation codes (C): FFTW, OMP
- Data visulatization codes (Python): Matplotlib, Plotly, Numpy, Pandas

## SCFT Calculations
Inputs
- Initial guess (seed): options for polymer field _w_(**r**) or polymer density distribution $\phi_p$(**r**)
- Polymer
  - Backbone charge fraction $\alpha$ (currently only homogeneous, though framework for multiblock has been built)
  - Flory--Huggins interaction parameter $\chi$ (currently only homogeneous, though framework for multiblock has been built)
  - Kuhn length $b$
  - Monomer volume $\nu$
- Brush grafting density $\sigma$
- Bulk salt concentration $c_\pm^b$
- System 
  - Box size $(L_x, L_y, L_z)$, discretization $(dx, dy, dz)$
  - Temperature, solution dielectric constant
- Convergence
  - Anderson Mixing number of histories to base next guess from $n_A$, mixing parameters $\sigma_1$ and $\sigma_2$

Outputs
- `ph*.dat`: polymer and solvent density distribution
- `el*.dat`: electrostatic potential, ion distributions, polymeric charge distribution
- `W*.dat`: polymer, solvent, and incompressibility field; electrostatic potential
- `printout.dat`: free energy and errors of the most recent iteration
- `it*.dat`: history of free energies and errors over calculation run

## Data visualization
Inputs: `ph*.dat` with corresponding box size, discretization, etc.

Outputs: Matplotlib, Plotly figures of polymer density heatmaps and isosurfaces


