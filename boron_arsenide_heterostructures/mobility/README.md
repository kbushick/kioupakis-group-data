# BAs Mobility
This subfolder contains the raw mobility data for strained and unstrained BAs calculated with EPW. The contents are as follows:
- converged\_strained contains the electron and hole mobility data computed on a 72x72x72 fine k and q grid at T = 100-500 K for the strained structure
- converged\_unstrained contains the electron and hole mobility data computed on a 72x72x72 fine k and q grid at T = 100-500 K for the unstrained structure
- strained\_convergence contains the electron and hole mobility data computed on a range of fine k and q grids at T = 300 K for the strained structure. The input and output epw2 files give an example for the 120x120x120 grid. The extrapolate\_mobilities.gnu file computes the mobility at a hypothetical infinitely fine mesh, which are reported in the paper. 
- unstrained\_convergence contains the electron and hole mobility data computed on a range of fine k and q grids at T = 300 K for the unstrained structure. The input and output epw2 files give an example for the 120x120x120 grid. The extrapolate\_mobilities.gnu file computes the mobility at a hypothetical infinitely fine mesh, which are reported in the paper.
- The liu\_temp data files contain data extracted from Phys. Rev. B 98, 081203 (2018), which are included in Figure 4 as a point of comparison. These data are for unstrained BAs and include spin orbit coupling. 
- temp\_plot.py generates Figure 4. 
- converge\_plot.py generates Figure S5.
