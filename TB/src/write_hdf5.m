nv = 2;
nc = 1;
h5create('myfile.h5','/crystal/alat',1)
h5write('myfile.h5','/crystal/alat',aXM1*sqrt(3))
h5create('myfile.h5','/crystal/avecs',size(unit_cell))
h5write('myfile.h5','/crystal/avecs',unit_cell)
h5create('myfile.h5','/crystal/centre',[1 3])
h5write('myfile.h5','/crystal/centre',(a1+a2)/2)
h5create('myfile.h5','/crystal/motif',[tot_natoms 3])
h5write('myfile.h5','/crystal/motif',[structure.x(:) structure.y(:) structure.z(:)])
h5create('myfile.h5','/crystal/n_atoms',1)
h5write('myfile.h5','/crystal/n_atoms',tot_natoms)
h5create('myfile.h5','/crystal/n_orbs',1)
h5write('myfile.h5','/crystal/n_orbs',tot_norbs/2)
h5create('myfile.h5','/crystal/orb_pattern',[1 tot_natoms])
ind1 = find(structure.name == "Mo");
ind2 = find(structure.name == "S");
orb_pattern(ind1) = 5;
orb_pattern(ind2) = 3;
h5write('myfile.h5','/crystal/orb_pattern',orb_pattern)
h5create('myfile.h5','/crystal/positions',size(mcell))
h5write('myfile.h5','/crystal/positions',mcell)
h5create('myfile.h5','/eigensystem/eigenvalues',[nv knum_tot])
h5write('myfile.h5','/eigensystem/eigenvalues',tb_bands(noccs-1:noccs,:))
h5create('myfile.h5','/eigensystem/eigenvectors_real',[nc tot_norbs knum_tot])
h5write('myfile.h5','/eigensystem/eigenvectors_real',real(tb_vecs(noccs+1,:,:)))
h5create('myfile.h5','/eigensystem/eigenvectors_imag',[nc tot_norbs knum_tot])
h5write('myfile.h5','/eigensystem/eigenvectors_imag',imag(tb_vecs(noccs+1,:,:)))
h5create('myfile.h5','/eigensystem/kgrid',[knum_tot 3])
h5write('myfile.h5','/eigensystem/kgrid',all_kpts)
h5create('myfile.h5','/eigensystem/n_k',1)
h5write('myfile.h5','/eigensystem/n_k',knum)
h5create('myfile.h5','/eigensystem/n_con',1)
h5write('myfile.h5','/eigensystem/n_con',nc)
h5create('myfile.h5','/eigensystem/n_val',1)
h5write('myfile.h5','/eigensystem/n_val',nv)