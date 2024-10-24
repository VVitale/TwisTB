function deformation = deformation_potential_porbs(deformation,coor,natoms,mcell,a0,nlayers)
def_avg = zeros(natoms,1);
ind_chalc = find(strcmp(coor.name,'S'));
if(isempty(ind_chalc))
    ind_chalc = find(strcmp(coor.name,'Se'));
end
natoms_chalc = length(ind_chalc);
for iiat = 1 : natoms_chalc
    iat = ind_chalc(iiat);
    iatom = [coor.x(iat) coor.y(iat) coor.z(iat)];
    nn = 0;
    nn_metal_atoms = zeros(3,3);
    layer = coor.layer(iat);
    % Area of hexagon composed of Chalcogen atoms
    S0 = sqrt(3)/4 * a0(layer)^2;
    h = (1.5 - 2*0.6226)*mcell(3,3);
    V0 = (S0 * h/2)/3;
    for jat = 1 : natoms
        if(layer == coor.layer(jat) && iat ~= jat)
            for ic = -1 : 1
                for jc = -1 : 1
                    jatom = [coor.x(jat) coor.y(jat) coor.z(jat)] + ic*mcell(1,:) + jc*mcell(2,:);
                    dist = iatom - jatom;
                    if(abs(norm(dist(1:2)) - a0(layer)/sqrt(3)) < 0.1)
                        nn = nn + 1;                            
                        if(nn > 3)
                            error('Too many neighbors')
                        end
                        nn_metal_atoms(nn,:) = dist;
                    end
                end
            end
        end
    end
    %nn_chalc_atoms(4,:) = nn_chalc_atoms(1,:);

    %rad = cart2pol(nn_chalc_atoms(:,1),nn_chalc_atoms(:,2));%poly2cw(nn_chalc_atoms(:,1),nn_chalc_atoms(:,2));
    %radWrapped = mod(rad,2*pi);
    %radWrapped(radWrapped==0 & rad>0) = 2*pi;
    %[~, sortIdx] = sort(radWrapped, 'descend');
    %sortedData = nn_chalc_atoms(sortIdx,:);  
    %S = polyarea(sortedData(:,1),sortedData(:,2));
    %deformation(iat) = ((abs((S_up-S0)) - abs(S_down-S0))/S0)*100; %tanh((dist-L/2)/4)*(1.0 -
    [K,V] = convhull([nn_metal_atoms;[0 0 0]]);
    deformation(iat) =(V-V0)/(V0); %tanh((dist-L/2)/4)*(1.0 -
end

%     chalc_coor = [coor.x(ind_chalc) coor.y(ind_chalc) coor.z(ind_chalc)];
%     supercell = [repmat(-mcell(1,:)-mcell(2,:),[natoms_chalc,1]);
%         repmat(-mcell(1,:),[natoms_chalc,1]);
%         repmat(-mcell(1,:)+mcell(2,:),[natoms_chalc,1]);
%         repmat(-mcell(2,:),[natoms_chalc,1]);
%         repmat([0 0 0],[natoms_chalc,1]);
%         repmat(mcell(2,:),[natoms_chalc,1]);
%         repmat(mcell(1,:)-mcell(2,:),[natoms_chalc,1]);
%         repmat(mcell(1,:),[natoms_chalc,1]);
%         repmat(mcell(1,:)+mcell(2,:),[natoms_chalc,1]);] + repmat(chalc_coor,[9,1]);
% index_atom_supercell = repmat([1:natoms_chalc]',[9,1]);
% average_deformation = average_radius(deformation(ind_chalc),natoms_chalc,chalc_coor,1.2*sqrt(3)*a0(1),supercell,index_atom_supercell);
% def_avg(ind_chalc) = average_deformation;
%figure
%scatter3(coor.x(ind_chalc),coor.y(ind_chalc),coor.z(ind_chalc),100,deformation(ind_chalc),'filled')
%colorbar

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
