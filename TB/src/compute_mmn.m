%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to compute the Mmn matrices
% from the TB eigenvectors in the 
% format required by Wannier90
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(mmn_fname,'w');
time = datestr(clock,'YYYY/mm/dd HH:MM:SS:FFF');
fprintf(fid,'%50s\n',join(['File generated with twisted bilayer studio. ',time]));
num_loc_bands = tot_norbs-num_exclude_bands;
fprintf(fid,'%i \t %i \t %i\n',num_loc_bands,knum_tot,num_nnkpt);
iik = 0;
loc_bands = setdiff([1:tot_norbs],exclude_bands);
% check the number of bands is correct
if(num_loc_bands ~= length(loc_bands))
	error('Wrong number of bands!')
end


% MAIN LOOP
for ikpt = 1 : knum_tot
	for innkp = 1 : num_nnkpt
		iik = iik + 1;
		% G vector connecting the k' vector in the first BZ to the actual periodic image
		gvec = nnkpt(iik,3)*mb1 + nnkpt(iik,4)*mb2 + nnkpt(iik,5)*mb3;
		% Compute e^(i G.r)
                eigr = exp(-1i.*(gvec*coords'));
		% index of k' point in the first BZ
		kb_index = nnkpt(iik,2);
		% Print 5 integers for each knum_tot*num_nnkpt block
		fprintf(fid,'%i \t %i \t %i \t %i \t %i\n',nnkpt(iik,:));
		% Print the real and imaginary parts of Mmn
                for mbnd = 1 : num_loc_bands
			for nbnd = 1 : num_loc_bands
				Mmn = dot(tb_vecs(:,loc_bands(nbnd),ikpt),tb_vecs(:,loc_bands(mbnd),kb_index).*eigr.');
				fprintf(fid,'\t %+8.12f   %+8.12f\n',real(Mmn),imag(Mmn));
			end
		end
	end
end
