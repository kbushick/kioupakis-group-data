# BAs Effective Mass Data
This subfolder contains the raw data used to calculate the effective masses of the carriers in BAs. The data (band\*.csv) is the same present in the gw\_so subfolder in the bands directory, recopied here for convinience. A helper bash script, `make_slices.sh` is provided to pare down these full bands (across the entire high symmetry path) into "slices" of data near the band extrema, necessary for calculating effective masses. These slices are organized in subfolders corresponding to the three different high symmetry directions of interest, G-X, G-K, and G-L. 

In the gamma\_\* subfolders, there are .gnu files that fit the band edges according to Equation 1 to obtain the effective masses. They also plot the band edges for visual verification. Note that masses are calculated and stored in the variable(s) m\[1-6\]. 

In the gamma\_x subfolder, there is also data for the transverse direction. slice9\_tran.csv contains the first 8 bands (with spin-orbit coupling) calculated along a direction perpindicular to G-X at the ∆ point. slice9\_tran\_limited.csv contains the same data, but only the x-step and lowest conduction band. Note that the x-units are different from this calculation, but this is accounted for in the find\_eff\_mass\_cbm\_transverse.gnu file.

Note that the different directions have different numbers of slices based on different degeneracies and the location of the conduction band minimum. 
