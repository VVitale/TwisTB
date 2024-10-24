%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to compute the Amn matrices
% from the TB eigenvectors in the
% format required by Wannier90
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(amn_fname,'w');
time = datestr(clock,'YYYY/mm/dd HH:MM:SS:FFF');
fprintf(fid,'%50s\n',join(['File generated with twisted bilayer studio. ',time]));
num_loc_bands = tot_norbs-num_exclude_bands;
if(~lsoc)
    fprintf(fid,'%i \t %i \t %i\n',num_loc_bands,knum_tot,num_proj);
else
    fprintf(fid,'%i \t %i \t %i\n',num_loc_bands,knum_tot,num_sproj);
end
iik = 0;
loc_bands = setdiff([1:tot_norbs],exclude_bands);

% check the number of bands is correct
if(num_loc_bands ~= length(loc_bands))
    error('Wrong number of bands!')
end

if(~lsoc)
    if(num_proj == tot_norbs)
        % Make sure we have exactly 5d orbitals per metal atom
        % and 3p orbitals per chalcogen atom
        if(nnz(ismember(proj_qnum(:,1),1))~=tot_norbs/11*6)
            error('number of p orbitals in projection block not correct!')
        end
        if(nnz(ismember(proj_qnum(:,1),2))~=tot_norbs/11*5)
            error('number of d orbitals in projection block not correct!')
        end
        % Assign each projection to an orbital in the TB basis
        proj2orb = zeros(num_proj,1);
        for iproj = 1 : num_proj
            for jorb = 1 : tot_norbs
                if(norm(proj_centres(iproj,:)-coords(jorb,:))<1e-4 &&...
                        proj_qnum(iproj,1) == orbitals(jorb).l && ...
                        proj_qnum(iproj,2) == orbitals(jorb).m)
                    proj2orb(iproj) = jorb;
                end
            end
        end
    end

    % MAIN LOOP
    if(num_proj == tot_norbs)
        if(num_loc_bands == num_proj)
            for ikpt = 1 : knum_tot
                for iwann = 1 : num_proj
                    for iband = 1 : num_proj
                        amn = conj(tb_vecs(proj2orb(iwann),iband,ikpt))*exp(-1.i*dot(all_kpts(ik,:),coords(proj2orb(iwann),:)));
                        fprintf(fid,'%i \t %i \t %i \t %+2.8f \t %+2.8f\n', iband, proj2orb(iwann), ikpt, real(amn), imag(amn));
                    end
                end
            end
        end
    elseif(num_proj < tot_norbs && num_proj > 0)
        % NB: only s orbitals are allowed in this case!
        % All these checks should be moved in w90_nnkp
        if(any(proj_qnum(:,1)~=0))
      		 error('Only s orbitals can be specified if num_proj not the same as total number of orbitals!')
        end
        rcut_thr = 1e-5;
        for ikpt = 1 : knum_tot
            for iband = 1 : num_loc_bands
                % Run over all orbitals and compute integrals, assuming each orbital to be
                % delta function
                for iwann = 1 : num_proj
                    % Compute cutoff radius
                    zona = proj_zona(iwann);
                    rcut = log(rcut_thr * sqrt(4*pi) / (2*zona))*(-1.0/zona)
                    if(rcut > min(norm(mcell(1,:)),norm(mcell(2,:)))/2)
                        error('Projection function too spread out. Increase zona value!')
                    end
                    % compute integrals including only those orbitals whose distance from the center is less then rcut
                    amn = complex(0.0,0.0);
                    for jorb = 1 : tot_norbs
                        for icell = -1 : 1
                            for jcell = -1 : 1
                                Rtauj = coords(jorb,:)+icell*mcell(1,:)+jcell*mcell(2,:);
                                if(norm(proj_centres(iwann,:)-Rtauj)<rcut)
                                    amn = amn + exp(1i*dot(all_kpts(ikpt,:),Rtauj))*tb_vecs(jorb,loc_bands(iwann),ikpt);
                                end
                            end
                        end
                    end
                    fprintf(fid,'%i \t %i \t %i \t %+2.8f \t %+2.8f\n', iband, iwann, ikpt, real(amn), imag(amn));
                end
            end
        end
    end
else
    if(num_sproj == tot_norbs)
        % Make sure we have exactly 5d orbitals per metal atom
        % and 3p orbitals per
        if(nnz(ismember(sproj_qnum(:,1),1))~=tot_norbs/11*6)
            error('number of p orbitals in projection block not correct!')
        end
        if(nnz(ismember(sproj_qnum(:,1),2))~=tot_norbs/11*5)
            error('number of d orbitals in projection block not correct!')
        end
        % Assign each projection to an orbital in the TB basis
        sproj2orb = zeros(num_sproj,1);
        for iproj = 1 : num_sproj
            for jorb = 1 : tot_norbs
                if(norm(sproj_centres(iproj,:)-coords(jorb,:))<1e-4 &&...
                        sproj_qnum(iproj,1) == orbitals(jorb).l && ...
                        sproj_qnum(iproj,2) == orbitals(jorb).m && ...
                        sproj_spin(iproj) == orbitals(jorb).Spin)
                    sproj2orb(iproj) = jorb;
                end
            end
        end
    end

    % MAIN LOOP
    if(num_sproj == tot_norbs)
        if(num_loc_bands == num_sproj)
            for ikpt = 1 : knum_tot
                for iwann = 1 : num_sproj
                    for iband = 1 : num_sproj
                        amn = conj(tb_vecs(sproj2orb(iwann),iband,ikpt))*exp(-1.i*dot(all_kpts(ik,:),coords(sproj2orb(iwann),:)));
                        fprintf(fid,'%i \t %i \t %i \t %+2.8f \t %+2.8f\n', iband, sproj2orb(iwann), ikpt, real(amn), imag(amn));
                    end
                end
            end
        end

    elseif(num_sproj < tot_norbs && num_sproj > 0)
        % NB: only s orbitals are allowed in this case!
        % All these checks should be moved in w90_nnkp
        if(any(sproj_qnum(:,1)~=0))
      		 error('Only s orbitals can be specified if num_proj not the same as total number of orbitals!')
        end
        rcut_thr = 1e-5;
        for ikpt = 1 : knum_tot
            for iband = 1 : num_loc_bands
                % Run over all orbitals and compute integrals, assuming each orbital to be
                % delta function
                for iwann = 1 : num_sproj
                    % Compute cutoff radius
                    zona = sproj_zona(iwann);
                    rcut = log(rcut_thr * sqrt(4*pi) / (2*zona))*(-1.0/zona);
                    if(rcut > min(norm(mcell(1,:)),norm(mcell(2,:)))/2)
                        error('Projection function too spread out. Increase zona value!')
                    end
                    % compute integrals including only those orbitals whose distance from the center is less then rcut
                    amn = complex(0.0,0.0);
                    for jorb = 1 : tot_norbs
                        for icell = -1 : 1
                            for jcell = -1 : 1
                                Rtauj = coords(jorb,:)+icell*mcell(1,:)+jcell*mcell(2,:);
                                if(norm(sproj_centres(iwann,:)-Rtauj)<rcut)% && (orbitals(jorb).Spin == sproj_spin(iwann)))
                                    amn = amn + exp(1i*dot(all_kpts(ikpt,:),Rtauj))*tb_vecs(jorb,loc_bands(iwann),ikpt);
                                end
                            end
                        end
                    end
                    fprintf(fid,'%i \t %i \t %i \t %+2.8f \t %+2.8f\n', iband, iwann, ikpt, real(amn), imag(amn));
                end
            end
        end
    end
end
