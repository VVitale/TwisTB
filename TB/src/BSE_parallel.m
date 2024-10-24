%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to construct the BSE                           %%
%% Hamiltonian and diagonalise it                          %%
%% in tight-binding basis given in :                       %%
%% (1) Ridolfi et al. PRB 97, 205409 (2018)                %%
%% and (2) MacDonald et al. PRB 91, 075310                 %%
%% (2015)                                                  %%
%% - INPUT VARIABLES                                       %%
%%   vn, cn : Number of valence and conduction bands       %%
%%   Cv, CC : Coeffients of TB for valence and cond bands  %%
%%   Ev, Ec : Eigenvalues of TB                            %%
%%   knum_tot : total number of k-points in k-grid         %%
%%   all_kpts : Cartesian coordinates of k-points          %%
%%   a : lattice constant                                  %%
%%   gradH_x : Gradient wrt k_x of Hamiltonian matrix      %%
%%   outfname : Name of output file (string)               %%
%%   mcell : Lattice cell vectors                          %%
%%   nat : Total number of atoms                           %%
%%   structure : Cartesian coordinates of all atoms        %%
%%   elho_int : Wheter to include el-ho interaction        %%
%%   dL : broadening of Lorentzian for plotting            %%
%%                                                         %%
%% - OUTPUT                                                %%
%%   BSE_eig, D1 : Eigenvalue of BSE w/o interaction       %%
%%   omega : array with h\omega/2\pi for plotting          %%
%%   Resigma : Real part of sigma_xx                       %%
%%   Acvk : Eigenvectors of BSE                            %%
%%   s : Oscillator strenght                               %%
%%   indc, indv, indk : Indices for c,v,k in BSE Hamilton. %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                         %%
%% Written by Valerio Vitale                               %%
%%                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [BSE_eig,omega,Resigma,Acvk,s,indc,indv,indk] = BSE_parallel2(vn,cn,Cv,Cc,Ev,Ec,knum_tot,all_kpts,...
                                                 a,Hmat,gradH_x,orb_xcoords,outfname,mcell,index,nat,...
                                                 structure,elho_int,dL)

% Dimension of BSE Hamiltonian 
dim = vn*cn*knum_tot;
% Number of BSE eigenvalues to compute (for the moment compute all of them)
neig = 50;% vn*cn*knum_tot;
% Volume of unit cell
V = sqrt(3)*a^2/2*knum_tot;
hbar = 6.582119569*10^(-16);
bse_thr = 1E-10;

% Fix gauge of Cv and Cc by imposing their sum to be
% a real number
for iv = 1 : vn
   for ik = 1 : knum_tot
      s = sum(Cv(:,iv,ik));
      phi = atan2(imag(s),real(s));
      phase = exp(-1i*phi);
      Cv(:,iv,ik) = phase.*Cv(:,iv,ik);
      if(imag(sum(Cv(:,iv,ik)))>10^(-10))
          [iv,ik]
          error('The sum of the valence band coefficients is not real.')
      end
   end
end
for ic = 1 : cn
   for ik = 1 : knum_tot
      s = sum(Cc(:,ic,ik));
      phi = atan2(imag(s),real(s));
      phase = exp(-1i*phi);
      Cc(:,ic,ik) = phase*Cc(:,ic,ik);
      if(imag(sum(Cc(:,ic,ik)))>10^(-10))
          [ic,ik]
          error('The sum of the conduction band coefficients is not real.')
      end
   end
end

% Define indices for conduction and valence bands and k-points
ind = zeros(dim,3);
idim = 0;
for ic = 1 : cn
   for iv = 1 : vn
       for ik = 1 : knum_tot
          idim = idim + 1;
          ind(idim,:) = [ic,iv,ik];
       end
   end
end

indc = ind(:,1);
indv = ind(:,2);
indk = ind(:,3);

omega = linspace(0.5,3.5,1000);
Resigma = zeros(length(omega),1);

% Compute absorbance spectrum with interaction
BSE_eig = zeros(neig,1);
if(elho_int)
   hold on
   clear D Resigma

   fprintf('--> Starting parallel construction of BSE Hamiltonian ...\n')

   % Compute lattice Fourier transform of Keldysh potential
   ukkp = complex(zeros(nat));
   u0 = complex(zeros(nat));
   u0 = screenC3([0,0,0],sqrt(knum_tot),nat,structure,mcell,a);

   % Compute diagonal part
   diag_H_BSE = zeros(dim,1);
   parfor idim = 1 : dim
      ic = indc(idim);
      iv = indv(idim);
      ik = indk(idim);
      Icc = zeros(nat,1);
      Ivv = zeros(nat,1);
      for inat = 1 : nat
              Icc(inat) = dot(Cc(index{inat},ic,ik),Cc(index{inat},ic,ik));
              Ivv(inat) = dot(Cv(index{inat},iv,ik),Cv(index{inat},iv,ik));
      end
      diag_H_BSE(idim) = Ec(ic,ik) - Ev(iv,ik) + Icc'*(u0*Ivv)/knum_tot;
   end

   shift = 2/3*mcell(1,:)+1/3*mcell(2,:);
   Ni = -(sqrt(knum_tot)-1)/2;
   Nf = -Ni;
   dN = 1;
   kn = 0;
   VR_ij = zeros(nat,nat,knum_tot);
   X_ij = zeros(nat,nat,knum_tot);
   Y_ij = zeros(nat,nat,knum_tot);
   Z_ij = zeros(nat,nat,knum_tot);
   R = zeros(knum_tot,3);
   imn = 0;
   for im = Ni : Nf
       for jn = Ni : Nf
           imn = imn + 1;
           R(imn,:) = im*mcell(1,:) + jn*mcell(2,:);
       end
   end
 
   parfor iR = 1 : knum_tot
       for inat = 1 : nat
           ti = [structure.x(inat),structure.y(inat),structure.z(inat)] - shift;
           for jnat = 1 : nat
               tj = [structure.x(jnat),structure.y(jnat),structure.z(jnat)] - shift;
               dist_ij = norm(R(iR,:) + tj - ti);
               VR_ij(inat,jnat,iR) = keldysh_pot(dist_ij,inat,jnat,a);
               X_ij(inat,jnat,iR) = R(iR,1) + tj(1) - ti(1);
               Y_ij(inat,jnat,iR) = R(iR,2) + tj(2) - ti(2);
               Z_ij(inat,jnat,iR) = R(iR,3) + tj(3) - ti(3);
           end
       end
   end

   %kx = repmat(all_kpts(:,1),[1, knum_tot]);
   %ky = repmat(all_kpts(:,2),[1, knum_tot]);
   %kz = repmat(all_kpts(:,3),[1, knum_tot]);

   %qx = kx - all_kpts(:,1)';
   %qy = ky - all_kpts(:,2)';
   %qz = kz - all_kpts(:,3)';

   %ukkp_loc = zeros(nat);
   ukkp = cell(knum_tot,knum_tot);
   for ik = 1 : knum_tot
      for ikp = 1 : knum_tot
         q = all_kpts(ik,:) - all_kpts(ikp,:);
         ukkp{ik,ikp} = zeros(nat);
         for iR = 1 : knum_tot
            qdotRij = q(1)*X_ij(:,:,iR) + q(2)*Y_ij(:,:,iR) + q(3)*Z_ij(:,:,iR);
            M = exp(1i*qdotRij);
            ukkp{ik,ikp} = ukkp{ik,ikp} +  M.*VR_ij(:,:,iR);
         end
      end
   end

   tic
    %spmd
    % Build BSE Hamiltonian upper triangular part
    %N1 = dim; N2 = dim;
    %globalSize_H_BSE = [N1 N2];
    %codistr_H_BSE = codistributor1d(2, codistributor1d.unsetPartition, globalSize_H_BSE);
    %localSize_H_BSE = [N1, codistr_H_BSE.Partition(labindex)];
    %H_BSE_loc = zeros(localSize_H_BSE);
    %H_BSE_loc_sparse = sparse(H_BSE_loc);
    %globalInd = codistr_H_BSE.globalIndices(2);

   H_BSE = complex(zeros(dim));
   %for jdim_loc = 1 : length(globalInd)
   %    jdim = globalInd(jdim_loc);
   parfor jdim = 1 : dim
       icp = indc(jdim);
       ivp = indv(jdim);
       ikp = indk(jdim);
       v = complex(zeros(dim,1));
       for idim = jdim + 1 : dim
          ic = indc(idim);
          iv = indv(idim);
          ik = indk(idim);

          Icc = zeros(nat,1);
          Ivv = zeros(nat,1);

          for inat = 1 : nat
                  Icc(inat) = dot(Cc(index{inat},icp,ikp),Cc(index{inat},ic,ik));
                  Ivv(inat) = dot(Cv(index{inat},ivp,ikp),Cv(index{inat},iv,ik));
          end
          v(idim)= Icc'*(ukkp{ik,ikp}*Ivv)/knum_tot;
      end
      v(abs(v) < bse_thr ) = complex(0.0,0.0);
      %H_BSE_loc_sparse(:,jdim) = sparse(v(:));
      H_BSE(:,jdim) = sparse(v(:));
      fprintf('Row %i completed \n', jdim) 
   end
   %H_BSE_global = codistributed.build(H_BSE_loc_sparse, codistr_H_BSE);
   %end
   toc

   %H_BSE = gather(H_BSE_global);
   H_BSE = H_BSE + H_BSE' + sparse(diag(diag_H_BSE));
   %H_BSE = H_BSE + sparse(diag(diag_H_BSE));
   fprintf('done.\n')
   clear Icc Ivv H_BSE_global
   
   % Check H_BSE is symmetric
   if(norm(full(H_BSE)' - full(H_BSE)) > 1E-8)
      %norm(H_BSE' - H_BSE)
      error('ERROR: hole-electron matrix is not Hermitian. Aborting!')
   end

   % Diagonalise BSE Hamiltonian
   fprintf('--> Starting diagonalization of BSE Hamiltonian ...\n')
   [Acvk,BSE_eig] = eigs(H_BSE,neig,'SR');
   BSE_eig = diag(BSE_eig);
   
   fprintf('done.\n')
   
   % In neig != dim pad
   if(dim-neig > 0)
	   padding = dim-neig;
   else
	   padding = neig-dim;
   end
   Acvk = [Acvk,eye(dim,padding)];
   BSE_eig = [BSE_eig;1000*ones(padding,1)];
   

   % Compute absorbance spectrum 
   s = zeros(neig,1);
   sumM = 0;
   for idim = 1 : neig
      sumM =0;
      for jdim = 1 : dim
          sumM = sumM + Acvk(jdim,idim)*(Cv(:,indv(jdim),indk(jdim))'*...
		  (gradH_x{indk(jdim)} + 1i*orb_xcoords*Hmat{indk(jdim)})*...
		  Cc(:,indc(jdim),indk(jdim)));
      end
      s(idim) = abs(sumM)^2;
   end
   for iom = 1 : length(omega)
      L = zeros(neig,1);
      for idim = 1 : neig
         L(idim) = lorentzian(omega(iom),BSE_eig(idim),dL);
      end
      Resigma(iom) = 4*pi/(V*omega(iom))*dot(s,L);
   end
   
   % Plot absorbance spectrum
   plot(omega,Resigma,'r')

   % Save absorbance spectrum
   saveas(gca,join([outfname,'Re_sigma_xx.png']))
   file_ID = fopen(join([outfname,'_OpticalCond_with_inter.dat']),'w');
   for iom = 1 : length(omega)
       fprintf(file_ID,'%2.8f %4.8f\n',omega(iom),Resigma(iom));
   end
   fclose(file_ID);
   
   % Compute Joint density of states
   E = linspace(0,6,1000);
   rho2 = zeros(size(E,2),1);
   for ie = 1 : length(rho2)
      for idim = 1 : dim
         %rho1(ie) = rho1(ie) + lorentzian(E(ie),D1(idim),dL);
         rho2(ie) = rho2(ie) + lorentzian(E(ie),BSE_eig(idim),dL);

      end
   end
   
   rho2 = rho2/cn/vn/knum_tot;
   % Plot Joint density of states
   figure
   plot(E,rho2,'r')
   saveas(gca,join([outfname,'_JDOS.png']))
   file_ID = fopen(join([outfname,'_JDOS_with_inter.dat']),'w');
   for iom = 1 : length(E)
       fprintf(file_ID,'%2.8f %4.8f\n',E(iom),rho2(iom));
   end

   fclose(file_ID);
   % delete(gcp('nocreate'))
   %clear Acvk Cc Cv 
   clear H_BSE rho E L rho1 rho2
   clear Ec Ev
end

end
