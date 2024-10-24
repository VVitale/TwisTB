function deformation = deformation_potential_dorbs(deformation,coor,natoms,mcell,a0,...
    strain,plot_surf)
def_avg = zeros(natoms,1);
ind_met = find(strcmp(coor.name,'Mo'));
if(isempty(ind_met))
    ind_met = find(strcmp(coor.name,'W'));
end
natoms_met = length(ind_met);
for iiat = 1 : natoms_met
    iat = ind_met(iiat);
    iatom = [coor.x(iat) coor.y(iat) coor.z(iat)];
    nn_up = 0;
    nn_down = 0;
    chalc_atoms_up = zeros(3,3);
    chalc_atoms_down = zeros(3,3);
    layer = coor.layer(iat);
    % Area of triangle composed of Chalcoge atoms
    S0 = sqrt(3)/4 * (a0(layer)/strain(layer))^2;
    h = (1.5 - 2*0.6226)*mcell(3,3);
    V0 = (S0 * h/2)/3;
    for jat = 1 : natoms
        if(layer == coor.layer(jat) && iat ~= jat)
            for ic = -1 : 1
                for jc = -1 : 1
                    jatom = [coor.x(jat) coor.y(jat) coor.z(jat)] + ic*mcell(1,:) + jc*mcell(2,:);
                    dist = iatom - jatom;
                    if(abs(norm(dist(1:2)) - a0(layer)./sqrt(3)) < 0.1)
                        %norm(dist(1:2))*sqrt(3)/(a0(layer)/strain(layer)) - 1
            	         if(dist(3) > 0)
                             nn_up = nn_up + 1;
                             if(nn_up > 3)
                                 error('Too many neighbors')
                             end
                             chalc_atoms_up(nn_up,:) = dist;
                	 elseif(dist(3) < 0)
                             nn_down = nn_down + 1;
                             if(nn_down > 3)
                                 error('Too many neighbors')
                             end
                             chalc_atoms_down(nn_down,:) = dist;
                	 end
                    end
                end
            end
        end
    end
    pts_up = [chalc_atoms_up;[0 0 0]];
    pts_down = [[0 0 0];chalc_atoms_down];
    [K_up,V_up] = convhull(pts_up(:,1),pts_up(:,2),pts_up(:,3));
    [K_down,V_down] = convhull(pts_down(:,1),pts_down(:,2),pts_down(:,3));
    %pts = [chalc_atoms_up;chalc_atoms_down];
    %[K,V] = convhull(pts(:,1),pts(:,2),pts(:,3));
    deformation(iat) = ((V_down-V0)+(V_up-V0))/(2*V0);
     if(exist("plot_surf",'var') ...
              && plot_surf && iiat==1)
             figure
             trisurf(K_down,pts_down(:,1),pts_down(:,2),pts_down(:,3),'Facecolor','cyan');
             hold on
             plot3(pts_down(:,1),pts_down(:,2),pts_down(:,3),'.r')
     end
end

%met_coor = [coor.x(ind_met) coor.y(ind_met) coor.z(ind_met)];
%supercell = [repmat(-mcell(1,:)-mcell(2,:),[natoms_met,1]);
%    repmat(-mcell(1,:),[natoms_met,1]);
%    repmat(-mcell(1,:)+mcell(2,:),[natoms_met,1]);
%    repmat(-mcell(2,:),[natoms_met,1]);
%    repmat([0 0 0],[natoms_met,1]);
%    repmat(mcell(2,:),[natoms_met,1]);
%    repmat(mcell(1,:)-mcell(2,:),[natoms_met,1]);
%    repmat(mcell(1,:),[natoms_met,1]);
%    repmat(mcell(1,:)+mcell(2,:),[natoms_met,1]);] + repmat(met_coor,[9,1]);
%index_atom_supercell = repmat([1:natoms_met]',[9,1]);
%average_deformation = average_radius(deformation(ind_met),natoms_met,met_coor,1.2*sqrt(3)*a0(1),supercell,index_atom_supercell);
%def_avg(ind_met) = average_deformation;
figure
scatter3(coor.x(ind_met),coor.y(ind_met),coor.z(ind_met),100,deformation(ind_met),'filled')
colorbar

% Computes macroscopic average over a lattice
    function value_on_lattice = average_radius(chern_op,num_atoms,at_pos,cutoff,at_pos_supercell,...
    		index_atom_supercell)
    	for iat = 1 : num_atoms
    		dist_all = sum((at_pos(iat,:) - at_pos_supercell).^2,2);
    		in_cutoff = find(dist_all <= cutoff^2);
    		value_on_lattice(iat) = sum(chern_op(index_atom_supercell(in_cutoff)))./length(in_cutoff);
    	end
    end
end
