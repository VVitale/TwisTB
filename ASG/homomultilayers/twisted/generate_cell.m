function [pos,k1,id] = generate_cell(alatL,ialatL,a,c,u,tot_nat,nlayers,theta,alat1,orientations,...
    stackings,atomlabels)
   % number of supercell for the search
   period = norm(alatL,1);
   tnat = floor(period/a);
   pos = -10^5*ones(tot_nat,3);
   size(pos)
   id(1,tot_nat) = "X";
   delta =  1*10^(-5);
   
%    disp(' ')
%    msg = ['Generating Cartesian atomic coordinates for moirè lattice with a twist angle of ',str_theta,' ...'];
%    disp(msg)
   
   k1=1;
   
   for ilayer = 1 : nlayers
       st=sin(theta(ilayer));
       ct=cos(theta(ilayer));
       Rtheta = [ct -st; st ct]
       for nj = -3*tnat:3*tnat
           for ni = -3*tnat:3*tnat
               r1 = ni*alat1(:,1) + nj*alat1(:,2);
               r2 = r1 + 1/3*alat1(:,1) + 1/3*alat1(:,2);
               %r3 = r2;
               f13 = 1/3;
               f23 = 2/3;
               if(strcmp(orientations{ilayer},'r'))
                   r1 = -r1 + 1/3*(alat1(:,1)+alat1(:,2));
                   r2 = -r2 + 1/3*(alat1(:,1)+alat1(:,2));
                   %r3 = -r2;
                   f13 = 2/3;
                   f23 = 1/3;
               end
               if(strcmp(stackings{ilayer},'t0'))
                   r1r = Rtheta*r1;
                   r2r = Rtheta*r2;
                   r3r = Rtheta*r2;
               elseif(strcmp(stackings{ilayer},'t13'))
                   r1r = Rtheta*(r1+f13*(alat1(:,1)+alat1(:,2)));
                   r2r = Rtheta*(r2+f13*(alat1(:,1)+alat1(:,2)));
                   r3r = r2r;
               elseif(strcmp(stackings{ilayer},'t23'))
                   r1r = Rtheta*(r1+f23*(alat1(:,1)+alat1(:,2)));
                   r2r = Rtheta*(r2+f23*(alat1(:,1)+alat1(:,2)));
                   r3r = r2r;
               end
               l1 = ialatL*r1r;
               l2 = ialatL*r2r;
               l3 = ialatL*r3r;
               if (l1(1) < (1-delta) && l1(1) >= -delta && l1(2) < (1-delta) && l1(2) >= -delta)
                   pos(k1,1)=r1r(1);
                   pos(k1,2)=r1r(2);
                   pos(k1,3)=((ilayer-1)*0.5 + 0.75)*c;
                   id(k1) = atomlabels{(ilayer-1)*3 + 1};
                   k1 = k1 + 1;
               end
               if (l2(1) < (1-delta) && l2(1) >= -delta && l2(2) < (1-delta) && l2(2) >= -delta)
                   pos(k1,1)=r2r(1);
                   pos(k1,2)=r2r(2);
                   pos(k1,3)= ((ilayer-1)*0.5 + u)*c;
                   id(k1) = atomlabels{(ilayer-1)*3 + 2};
                   k1 = k1 + 1;
              % end
              % if (l3(1) < (1-delta) && l3(1) >= -delta && l3(2) < (1-delta) && l3(2) >= -delta)
                   pos(k1,1)=r3r(1);
                   pos(k1,2)=r3r(2);
                   pos(k1,3)=((ilayer-1)*0.5+1.5-u)*c;
                   id(k1) = atomlabels{(ilayer-1)*3 + 3};
                   k1 = k1 + 1;
               end
           end
       end
   end
   

end