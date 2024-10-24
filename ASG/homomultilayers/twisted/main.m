 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to find a set of commensurate moire lattices 
% and the atomic coordinates for twisted homobilayers.
% The twisted bilayer is assumed to be made of 
% a material whose conventional cell is hexagonal, 
% e.g. graphene, hexagonal boron nitride (hBN) and 
% Transition-metal dichalcogenides (TMDs)
%
%
% INPUT: 
%  n      - Integer parameter that defines Moire cell,
%           m=n+1 by default, see Ref. 1
%  a      - lattice parameter of 2D material
%  c      - Interlayer distance when starting 2D material
%           consists of a monolayer, e.g. graphene, OR 
%           disance between central atoms of 2D material 
%           for trilayers, e.g. distance 
%           between transition metal atoms in the two
%           TMD layers.
%  comb   - Define which combination of layers to consider
%           Possible choices are: 'graphene', 'hBN', 
%           'MoS2', 'WSe2', 'MoSe2', 'WS2'
% plot_f  - Format for plotting ('xyz' or 'xsf')
% kpath   - Path in k-space. Implemented paths are:
%           * K-Gamma-M-K' -> 'KMGKp'
%           * M-Gamma-K-M  -> 'MGKM'
%           * Gamma-M-K-Gamma -> 'GMKG'
%           * Gamma-K-M    -> 'GKM'
%           * K-M          -> 'KM'
%           * Monkhorst-Pack -> 'uniform'
% init_nk - Total number of kpoints for uniform MP grids OR 
%           number of k-points in the first segment of the 
%           chosen path
%
% OUTPUT:
%  <rootname>.xyz(xsf)  - Cartesian coordinates in XYZ 
%                         (XSF) format
%  positions.<rootname>.dat  - Cartesian coordinates & 
%                              Moire lattice vectors
%                              to be used as input file
%                              for TB_free.x
%  lammps_positions.<rootname>.dat - Lammps geometry
%                                    file. 
%                                    N.B. The Moire 
%                                    cell is rotated
%                                    such as one of 
%                                    the vector is
%                                    parallel to
%                                    on the x-axis
%  potential.<rootname>.dat - Potential file to be
%                             used as input in 
%                             TB_free.x. Only for hBN
%                             at the moment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Written by Valerio Vitale, October 2019
%
% Refs. 
%       1) Nano Lett. 10, 3, 804-808 (2010)
%       2) PRB 86, 155449 (2012) 
%
% version 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 11-band TB model PRB 92 205108 (2015)
% All parameters are in Ang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        a      d(M-X)   d(X-X)    c (sqrt(d(M-X)^2-1/3*a^2)/(3/4-u))    u
% MoS2  3.1824  2.41      3.13    12.29
% WS2   3.1817  2.42      3.14    12.32
% MoSe2 3.3174  2.54      3.34    12.90
% WSe2  3.3155  2.55      3.35    12.96

% Set from SW classical potential LAMMPs
% All parameters are in Ang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        a      c       u 
% MoS2  3.110   12.08  0.621
% WS2   3.129   12.12  0.622
% MoSe2 3.3096  12.84  0.621
% WSe2  3.289   12.74  0.621
    
% Parameters from our DFT calculations
% All parameters are in Ang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          a      c         u
% MoS2  3.183  12.0264    0.621
% WS2   3.182  12.1634    0.621
% MoSe2 3.315  12.9567    0.621
% WSe2  3.3145 13.0089    0.621

dft_OPTB88_qe_geom = [
 3.183  12.0264 
 3.315  12.9567 
 3.182  12.1634 
 3.3145 13.0089];

dft_PBE_D3_siesta_geom = [
 3.183  12.0941 
 3.315  12.6589 
 3.182  11.7851 
 3.3145 12.5227];

lammps_sw_geom = [
 3.117   12.00 
 3.311   12.84 
 3.127   12.21 
 3.288   12.82]; 

geom_data = dft_OPTB88_qe_geom; %lammps_sw_geom;
mat_params = cell(6,5);
%MoS2 imat = 1
mat_params{1,1} = 'Mo';mat_params{1,2} = 'S';mat_params{1,3} = geom_data(1,1);mat_params{1,4} = geom_data(1,2);mat_params{1,5} = 'TMD'; 
%MoSe2 imat = 2
mat_params{2,1} = 'Mo';mat_params{2,2} = 'Se';mat_params{2,3} = geom_data(2,1); mat_params{2,4} = geom_data(2,2);mat_params{2,5} = 'TMD';
%WS2 imat = 3
mat_params{3,1} = 'W';mat_params{3,2} = 'S';mat_params{3,3} = geom_data(3,1); mat_params{3,4} = geom_data(3,2);mat_params{3,5} = 'TMD';
%WSe2 imat = 4
mat_params{4,1} = 'W';mat_params{4,2} = 'Se';mat_params{4,3} = geom_data(4,1); mat_params{4,4} = geom_data(4,2);mat_params{4,5} = 'TMD';
%graphene imat = 5
mat_params{5,1} = 'C';mat_params{5,2} = 'C';mat_params{5,3} = 2.46; mat_params{5,4} = 3.30;mat_params{5,5} = 'graphene';
%hBN imat = 6
mat_params{6,1} = 'B';mat_params{5,2} = 'N';mat_params{5,3} = 2.504; mat_params{5,4} = 3.35/2;mat_params{6,5} = 'hBN'


%%%%%%%%%%%% BEGIN INPUT %%%%%%%%%%%%%
% Integer for constructing Moire supercell
n = 1;
% Material index (see table at the top)
imat = 2;
% Lattice parameter of monolayer
a = mat_params{imat,3};
%a = 3.182; %
% Interlayer distance
c = mat_params{imat,4};
%c = 12.1634/2;
% Info for writing Cartesian coordinates and potential to file
write_positions = true;
write_potential = false;
write_cart      = true;
% Whether to XYZ or XSF format
plot_f = 'xsf';
write_lammps_input = true;
% Initial phase
nlayers = 1;
orientations = {'u','r','r'};%,'r'};
translations = {'t0','t0','t0'};%,'t13'};
atomlabels = {mat_params{imat+2,1},join([mat_params{imat+2,2},'1']),join([mat_params{imat+2,2},'2']), ...
       mat_params{imat+2,1},join([mat_params{imat+2,2},'1']),join([mat_params{imat+2,2},'2']),...
       mat_params{imat,1},join([mat_params{imat,2},'1']),join([mat_params{imat,2},'2'])}  
% Chemical symbol for system
comb = mat_params{imat,5};
theta_vec = [0,0,1];%,1,1];

% Info for writing k-points and path in k-space
write_kpath = true;
kpath = 'GMKG';
init_nk = 10;
%%%%%%%%%%%% END INPUT %%%%%%%%%%%%%

switch comb
    case 'graphene'
        disp('Structure and input files for twisted bilayer graphene')
        disp(' ')
        [num,pos,pos2] = TBLG_cell(n,a,c,plot_f,kpath,init_nk, write_positions, write_potential, write_cart, write_lammps_input, write_kpath,phase);
    case 'hBN'
        disp('Structure and input files for twisted bilayer hBN')
        disp(' ')
        [num,pos,pos2] = TBLhBN_cell(n,a,c,plot_f,kpath,init_nk, write_positions, write_potential, write_cart, write_lammps_input, write_kpath,phase);
    case 'TMD'
        disp('Structure and input files for twisted bilayer TMD')
        disp(' ')
        [num,pos] = TMLTMD_cell(theta_vec,n,a,c,atomlabels,...
            plot_f,kpath,init_nk, write_positions, ...
            write_cart, write_lammps_input, write_kpath,translations,...
            orientations,nlayers);
    otherwise
        error('Material not supported. Available choices are: graphene, hBN, MoS2 and WSe2')
end
%end
disp('Done! Have a nice day :-)')
