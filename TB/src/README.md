To run the MATLAB code you will need an input file (default is 'input.dat').
Keywords to specify in the input file are described below

```
- task                   (I) 1: bandstructure  2: Chern number and Berry curvature   3: DOS  4: charge density   5: exciton spectra and JDOS
6: Plot wavefunctions at Gamma, M and K   7: Build Hamiltonian on uniform grid and write in hdf5 format   8: Compute eigenvectors and matrices for Wannier90
- restart                (L) Whether this is a restart calculation
- read_ham               (L) Whether to read Hamiltonian 
- geomfname              (S) Name of geometry file
- gw_par                 (L) Whether to use GW parameters
- d12                    (R) maximum distance between chalcogen atom on bottom layer and chalcogen atom on top layer for computing 
                              interlayer hopping (aXX2 = 5.0)
- nlayer                 (I) Number of layers
- tmdc                   (I) TMD material code --- 1: MoS2   2: MoSe2   3: WS2   4:WSe2, e.g. MoS2 monolayer tmdc: = 1, MoS2 bilayer: tmdc = 11, 
                              MoS2 trilayer tmdc : 111, WS2/MoS2 bilayer: tmdc = 13
- convention             (I) Convention for the unit cell vectors of the monolayer (a=sqrt(3)*aXM)
                                                          a1                   a2
                              convention 1 (120 deg) | (1, 0)*a   |  (-1/2,sqrt(3)/2)*a
                              convention 2 (60 deg)  | (1, 0)*a   |  (1/2,sqrt(3)/2)*a
- knum                   (I) Number of k-points in the first segment (Gamma-M) when task = 1. Number of k-points along one reciprocal
                              vector when task = {2,3,4,5,7}, in this case the total number of k-points is knum^2
- miniBZ                 (L) Whether to use Gamma-M-K-Gamma of the miniBZ of Gamma-M-K-Gamma of the BZ of the monolayer folded onto the miniBZ
- interlayer_int         (L) Whether to compute the interlayer interaction
- spin_orbit             (L) Whether to compute spin-orbit 
- reduced_workspace      (L) Whether to use eig or eigs in the diagonalisation
- num_eigs               (I) Number of conduction bands when reduced_workspace = true 
- eigvecs                (L) Whether to compute the eigenvectors
- save_workspace         (L) Whether to save the orbital class for restart calculations
- num_workers            (I,S) Number of workers in parallel diagonalisation or 'max' maximum number of workers allowed by local profile
- dos_spacing            (R) Energy spacing in eV for DOS plot (e.g. 0.005)
- gradient               (L) Whether to compute the x component of the Hamiltonian gradient
- Interpd                (L) Whether to compute the pz-dz2 hopping
- lambda                 (R) Shift of eigenspectrum, H = H + lambda*I
- bse_num_vbnd           (I) Number of valence bands in BSE
- bse_num_cbnd           (I) Number of conduction bands in BSE
- bse_int                (L) Whether to compute spectra with BSE interaction
- bse_broadening         (R) Broadening of the Lorentzian used to approximate the delta function in optical spectra
- bse_eig_plot           (I) Number of excitonic eigenvfunction to plot
- bse_irh                (I) Index of hole atomic position in the (moire) unit cell
- bse_serial             (L) Whether to run the BSE solver serially
- bse_shifted            (L) Whether to use a randomly shifted from (0,0) uniform k-grid for BSE
- write_ham              (L) Write Hamiltonian to hdf5 format
- ham_fname              (S) Name of hamiltonian hdf5 file
- flipped                (S) Whether a layer has been flipped to form a 2H stacking
- read_kpts              (L) Whether to read k-points
- ef_strength            (R) Electric-field strength in V/Ang
- onsite_moire           (R) Strength of onsite term proportional to local angle
- sixth_nn               (L) Whether to use the Hamiltonian with up to 6th nearest neighbour hoppings
- g1                     (R) Screened deformation potential for d orbitals
- w90_root               (S) root filename for Wannier90 files

R = real number
I = integer number
L = boolean (true, false)
S = string
```
