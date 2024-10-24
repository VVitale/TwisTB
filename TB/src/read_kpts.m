% Check kpts.in file exists in the path
if(exist('kpts.in','file'))
   data = importdata('kpts.in',' ')
   recL=[mb1;mb2;mb3];
   all_kpts = data*recL;
   knum_tot = size(all_kpts,1)
   scale_axis = zeros(knum_tot,1);
   for ik = 2 : knum_tot
       dk = norm(all_kpts(ik,:)-all_kpts(ik-1,:));
       scale_axis(ik) = scale_axis(ik-1) + dk;
   end
else
  error('Error main.m cannot find kpts.in. Aborting ...')
end
