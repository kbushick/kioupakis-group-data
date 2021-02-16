# BAs Optical Data, excitonic effects 
This subfolder contains the raw data used to compute the optical properties of BAs including excitonic effects, but without including indirect processes. 

- absorption\_eh.dat is the output from the absorption code in BerkeleyGW. The first column is energy, the second column is epsilon2, the third column is epsilon1, and the fourth column is the normalized DOS. 
- we still recalculate epsilon1 using Kramers-Kronig for consistency, using the fortran code in kk.F90
- calculate\_optical\_data.sh takes the epsilon2 vs energy data and performs a number of steps to calculate other important optical properties. The combination of epsilon1 (from KK) and epsilon2 are used to calculate n, k, and alpha, with the output saved in optical\_data\_eh.{dat,csv}. This is the final data used in the paper (Figure 3 and Figure 4).
