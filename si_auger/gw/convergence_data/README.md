# Si GW convergence
This subfolder contains the output data used to determine the converged parameters for GW of bulk Si using BerkeleyGW. We test the 4 parameters independently, keeping the other 3 at their "infinity" value while varying the fourth. For # epsilon bands we use 500, 601, 700, 800, and 999 as our values. For # sigma bands we use 512, 582, 648, 798, and 998 as our values (restricted by degeneracy). For the epsilon and sigma cutoff energies we use 20, 25, 30, 35, and 37 Ry.  

- In each folder, we provide the sigma\_hp.log file from these respective runs, which is renamed for the parameter that was altered. 
- These files are read by the provided python script, which is provided in BGW tutorial, and calculates differences between different k points, highlighting a few high symmetry points. 
- The output of the python script is the included graph saved as a png

