# r-GeO2 Data
This folder contains the raw data present in "Electron and hole mobility of rutile GeO2 from first principles: An ultrawide-bandgap semiconductor for power electronics" (Appl. Phys. Lett. 117, 182104 (2020); https://doi.org/10.1063/5.0033284). The data is organized into subfolders, and the jupyter notebook in this folder generates all of the figures from the paper. 

For EPW runs, one epw1 calculation (see intitial\_inputs) was run, and then the output from this run: vmedata.fmt, ksdata.fmt, epwdata.fmt, crystal.fmt, k\*map, \*epmatwp, \*ukk, and \*bvec, were used for different secondary calculations (mobility, elph coupling, self-energy). 

- The mobility subfolder corresponds to the mobility data presented in Figure 1.
- The electron\_phonon\_coupling subfolder corresponds to the elph ME data presented in Figure 2.
- The self\_energies subfolder corresponds to the self-energy data presented in Figure 3. 
- The initial\_inputs folder contains the pseudopotentials and input files used in multiple calculations for this work. 
- The external\_data contains data extracted from plots to compare GaN mobility values in the SI. 

The arxiv version of the paper can be found [here](https://arxiv.org/abs/1911.09750v2). 
