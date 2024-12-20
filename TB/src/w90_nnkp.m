%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Routine to read the .nnkp file 
% generated by Wannier90 and computes
% kpoints coordinates as 
% defined in the .win file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Written by Valerio Vitale
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Open .nnkp file
if(~exist(nnkp_fname,'file'))
	error('Reading w90.nnkp file does not exist!');
end
input = textread(nnkp_fname,'%s','delimiter','\n');

nlines = size(input,1);
i=0;found=false;
spinor_found=false;
proj_found=false;
for iline = 2 : nlines %skip first line
	if(~isempty(input{iline}))
		tline = textscan(input{iline},'%s%s','delimiter',':');
                tkey = deblank(tline{1});
		if(strcmp(tkey{1},'calc_only_A'))
			calc_only_A = strcmp(deblank(char(tline{2})),'T');
		elseif(strcmp(tkey,'begin real_lattice'))
			% initialise array
			real_latt = zeros(3,3);
			% check that the block is closed
			tkey2 = textscan(input{iline+4},'%s','delimiter','\n');
			if(~strcmp(tkey2{1},'end real_lattice'))
				error('Error in reading the real_lattice block. Block is not closed!')
			end
			for jline = 1 : 3
				lvec = textscan(input{iline+jline},'%f%f%f','delimiter','\t');
				real_latt(jline,:) = [lvec{1} lvec{2} lvec{3}];
			end
		elseif(strcmp(tkey,'begin recip_lattice'))
			% initialise array
			recip_latt = zeros(3,3);
			% check that the block is closed
			tkey2 = textscan(input{iline+4},'%s','delimiter','\n');
			if(~strcmp(tkey2{1},'end recip_lattice'))
				error('Error in reading the recip_lattice block. Block is not closed!')
			end
			for jline = 1 : 3
				lvec = textscan(input{iline+jline},'%f%f%f','delimiter','\t');
				recip_latt(jline,:) = [lvec{1} lvec{2} lvec{3}];
			end
		elseif(strcmp(tkey,'begin kpoints'))
        	        tkey2 = textscan(input{iline+1},'%f','delimiter','\n');
        	        num_kpts = tkey2{1};
			% initialise array
			frac_kpoints = zeros(num_kpts,3);
			% check that the block is closed
			tkey2 = textscan(input{iline+num_kpts+2},'%s','delimiter','\n');
			if(~strcmp(tkey2{1},'end kpoints'))
				error('Error in reading the kpoints block. Block is not closed!')
			end
			for ikpt = 1 : num_kpts
				lvec = textscan(input{iline+1+ikpt},'%f%f%f','delimiter','\t');
				frac_kpoints(ikpt,:) = [lvec{1} lvec{2} lvec{3}];
			end
        elseif(strcmp(tkey,'begin projections'))
            proj_found = true;
	    if(lsoc)
		    error('Found projections block but spin_orbit set to true!')
	    end
            tkey2 = textscan(input{iline+1},'%f','delimiter','\n');
            num_proj = tkey2{1};
		% initialise array
		proj_centres = zeros(num_proj,3);
            proj_qnum = zeros(num_proj,2);
            proj_r = zeros(num_proj,1);
            proj_zaxis = zeros(num_proj,3);
            proj_xaxis = zeros(num_proj,3);
            proj_zona = zeros(num_proj,1);
    	    % check that the block is closed
	    tkey2 = textscan(input{iline+num_proj*2+2},'%s','delimiter','\n');
	    if(~strcmp(tkey2{1},'end projections'))
		    error('Error in reading the projections block. Block is not closed!')
            end
            iiproj = 0;
	    for iproj = 1 : 2 : 2*num_proj
                iiproj = iiproj+1;
		lvec = textscan(input{iline+1+iproj},'%f%f%f%f%f%f','delimiter','\t');
		proj_centres(iiproj,:) = [lvec{1} lvec{2} lvec{3}]*real_latt;
                proj_qnum(iiproj,:) = [lvec{4} lvec{5}];
                proj_r(iiproj) = lvec{6};
                lvec = textscan(input{iline+2+iproj},'%f%f%f%f%f%f%f','delimiter','\t');
                proj_xaxis(iiproj,:) = [lvec{1} lvec{2} lvec{3}];
                proj_zaxis(iiproj,:) = [lvec{4} lvec{5} lvec{6}];
                proj_zona(iiproj) = lvec{7};
            end
         elseif(strcmp(tkey,'begin spinor_projections'))
            spinor_found = true;
            if(proj_found)
                error('Both spinor_projections and normal projections present!')
            end
	    if(~lsoc)
		    error('spinor_projections found but spin_orbit is false!')
	    end
            tkey2 = textscan(input{iline+1},'%f','delimiter','\n');
            num_sproj = tkey2{1};
	    % initialise array
	    sproj_centres = zeros(num_sproj,3);
            sproj_qnum = zeros(num_sproj,2);
            sproj_r = zeros(num_sproj,1);
            sproj_zaxis = zeros(num_sproj,3);
            sproj_xaxis = zeros(num_sproj,3);
            sproj_zona = zeros(num_sproj,1);
            sproj_spin = zeros(num_sproj,1);
            sproj_qntax = zeros(num_sproj,3);
	    % check that the block is closed
	    tkey2 = textscan(input{iline+num_sproj*3+2},'%s','delimiter','\n');
	    if(~strcmp(tkey2{1},'end spinor_projections'))
		    error('Error in reading the spinor_projections block. Block is not closed!')
	    end
            iisproj = 0;
	    for isproj = 1 : 3 : 3*num_sproj
    		    iisproj = iisproj+1;
		    lvec = textscan(input{iline+1+isproj},'%f%f%f%f%f%f','delimiter','\t');
		    sproj_centres(iisproj,:) = [lvec{1} lvec{2} lvec{3}]*real_latt;
    		    sproj_qnum(iisproj,:) = [lvec{4} lvec{5}];
		    sproj_r(iisproj) = lvec{6};
	    	    lvec = textscan(input{iline+2+isproj},'%f%f%f%f%f%f%f','delimiter','\t');
		    sproj_xaxis(iisproj,:) = [lvec{1} lvec{2} lvec{3}];
		    sproj_zaxis(iisproj,:) = [lvec{4} lvec{5} lvec{6}];
		    sproj_zona(iisproj) = lvec{7};                
		    lvec = textscan(input{iline+3+isproj},'%f%f%f%f','delimiter','\t');
		    sproj_spin(iisproj) = lvec{1};
		    sproj_qntax(iisproj,:) = [lvec{2} lvec{3} lvec{4}];
	    end   
        elseif(strcmp(tkey,'begin nnkpts'))
            tkey2 = textscan(input{iline+1},'%f','delimiter','\n');
            num_nnkpt = tkey2{1};
            % initialise array
            nnkpt = zeros(num_kpts*num_nnkpt,5);
            % check that the block is closed
            tkey2 = textscan(input{iline+num_kpts*num_nnkpt+2},'%s','delimiter','\n');
            if(~strcmp(tkey2{1},'end nnkpts'))
                error('Error in reading the nnkpt block. Block is not closed!')
            end
            iik = 0;
            for ikpt = 1 : num_kpts
                for innkpt = 1 : num_nnkpt
                    iik = iik + 1;
                    lvec = textscan(input{iline+1+iik},'%f%f%f%f%f%f','delimiter','\t');
                    nnkpt(iik,:) = [lvec{1} lvec{2} lvec{3} lvec{4} lvec{5}];
                end
            end
        elseif(strcmp(tkey,'begin exclude_bands'))
            tkey2 = textscan(input{iline+1},'%f','delimiter','\n');
            num_exclude_bands = tkey2{1};
            % initialise array
            exclude_bands = zeros(num_exclude_bands,1);
            % check that the block is closed
            tkey2 = textscan(input{iline+num_exclude_bands+2},'%s','delimiter','\n');
            if(~strcmp(tkey2{1},'end exclude_bands'))
                error('Error in reading the exclude_bands block. Block is not closed!')
            end
            for iexbnd = 1 : num_exclude_bands
                lvec = textscan(input{iline+1+iexbnd},'%f','delimiter','\n');
                exclude_bands(iexbnd) = lvec{1};         
            end
        end
    end
end

% Run checks
% Real Lattice
if(abs(dot(ma1,real_latt(1,:))/(norm(ma1)*norm(real_latt(1,:)))-1.0) >= 1e-8)
	error('Error in real lattice vector 1. Check the nnkp file!')             
end                                                                               
if(abs(dot(ma2,real_latt(2,:))/(norm(ma2)*norm(real_latt(2,:)))-1.0) >= 1e-8)
	error('Error in real lattice vector 2. Check the nnkp file!')             
end                                                                               
if(abs(dot(ma3,real_latt(3,:))/(norm(ma3)*norm(real_latt(3,:)))-1.0) >= 1e-8)
	warning('Error in real lattice vector 3. Check the nnkp file!')
end


% Reciprocal Lattice
% Reciprocal lattice vectors of moire cell
% These need to be recomputed here to avoid issues with different conventions
mv=abs(dot(ma1,cross(ma2,ma3)));
mb1=2*pi*cross(ma2,ma3)/mv;
mb2=2*pi*cross(ma3,ma1)/mv;
mb3=2*pi*cross(ma1,ma2)/mv;

if(abs(dot(mb1,recip_latt(1,:))/(norm(mb1)*norm(recip_latt(1,:)))-1.0) >= 1e-8)
	error('Error in reciprocal lattice vector 1. Check the nnkp fi le!')
end                                                                    
if(abs(dot(mb2,recip_latt(2,:))/(norm(mb2)*norm(recip_latt(2,:)))-1.0) >= 1e-8)
	error('Error in reciprocal lattice vector 2. Check the nnkp fi le!')
end                                                                    
if(abs(dot(mb3,recip_latt(3,:))/(norm(mb3)*norm(recip_latt(3,:)))-1.0) >= 1e-8)
	warning('Error in reciprocal lattice vector 3. Check the nnkp file!')
end

clear real_latt recip_latt

% Compute cartesian coordinates of reciprocal lattice vectors
recL = [mb1;mb2;mb3];
all_kpts = frac_kpoints*recL;
knum_tot = num_kpts;
scale_axis = zeros(knum_tot,1);
for ik = 2 : knum_tot
    dk = norm(all_kpts(ik,:)-all_kpts(ik-1,:));
    scale_axis(ik) = scale_axis(ik-1) + dk;
end


