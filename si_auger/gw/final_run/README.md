# Si GW 
This subfolder contains the raw data used to compute the quasiparticle corrections to the band structure of bulk Si using BerkeleyGW. Output files from the run are also included for reference, though the large data files used during the run (WFN, RHO, epsmat.h5, etc.) are not included.  

- The scf input and pseduo-potential are used to generate the charge density for the subsequent bands (nscf) runs in the wfn, wfnq, and wfn\_co steps.
- The wfn, wfnq, and wfn\_co folders all contain inputs to generate the different wavefunction and density files used in the epsilon and sigma steps. WFN and WFNq are used for the epsilon step, while WFN\_co (renamed WFN\_inner) and RHO are used in the sigma step. 
- epsilon contains the input for the epsilon step, as well as generated output from our run. 
- sigma contains input for both the sigma step within the BGW workflow. log, output, and dat files are also provided from the run. We also provide input for the post-processing sig2wan step, which generates the Si\_gw.eig file that is read by wannier90.
