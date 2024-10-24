%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to compute the .eig file
% from the TB eigenvalues in the 
% format required by Wannier90
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(eig_fname,'w');
time = datestr(clock,'YYYY/mm/dd HH:MM:SS:FFF');
num_loc_bands = tot_norbs-num_exclude_bands;
loc_bands = setdiff([1:tot_norbs],exclude_bands);
% check the number of bands is correct
if(num_loc_bands ~= length(loc_bands))
    error('Wrong number of bands!')
end

for ikpt = 1 : knum_tot
        for iband = 1 : num_loc_bands
		fprintf(fid,'%i \t %i \t %+2.8f\n',iband,ikpt,tb_bands(loc_bands(iband),ikpt));
	end
end
