# BAs Band Structure Data
This subfolder contains the raw data used to plot the band structure of BAs. The full lists of points and energies along high symmetry paths are available are present in the files starting with pz. Four sets of data are provided, corresponding to the four band structures presented in the paper. Two are computed with LDA, while two are computed using GW. Of each of these pairs, one does not include spin-orbit coupling and one does (denotes "so"). 

The python script used to generate the band structure reads each band separately, so a helper bash file `split_bands.sh` is also provided to separate out the full data files into individual band files. These are moved to their corresponding subdirectories: lda, lda_so, gw, and gw_so for clarity. 

Note that the data is plotted along the high-symmetry path: Gamma-X-W-K-Gamma-L-U-W-L-K/U-X. These high-symmetry points have x-coordinates: 0.00000, 1.31513, 1.97270, 2.43767, 3.83258, 4.97152, 5.77687, 6.24184, 7.17178, 7.97714, and 8.44211, respectively. 
