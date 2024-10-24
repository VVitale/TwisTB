function hopping = assign_hoppings_6thnn_strain(a1,a2,delta,neigh_shell,layer,...
    sym_sign,ir,jr,type1,type2,type3,type4,type5,type7,type8,type10,type12,...
    type14,type16,typem16,type18,typem18,type19,type20,typem21)%,del1,del2,...
    %del3,del4,del5,del6,del7,del8,del9,del10,del11,del12,del13,del14,del15,del16,...
    %del17,del18,del19,del20,del21)
hopping = zeros(size(delta,1),1);
tol = 0.3;
Lambda = zeros(11);
lu=0;lv=0;
for ind1 = 1 : 11
    if(ind1 == 1 || ind1 == 2 || ind1 == 6 || ind1 == 7 || ind1 == 8)
    	lu = 2;
    else
    	lu = 1;
    end
    for ind2 = 1 : 11
    	if(ind2 == 1 || ind2 == 2 || ind2 == 6 || ind2 == 7 || ind2 == 8)
            lv = 2;
        else
            lv = 1;
    	end
    	Lambda(ind1,ind2) = lu + lv + 1;
    end
end
Lambda = Lambda./5;

drr0 = 0.0;
if (neigh_shell == 1)
    del4=-(2*a1+a2)/3;
    del5=(a1+2*a2)/3;
    del6=(a1-a2)/3;
    del4_len = norm(del4);
    del5_len = norm(del5);
    del6_len = norm(del6);
    for ih = 1 : length(hopping)
	delta_len = norm(delta(ih,:));
        d4=dot(delta(ih,:),del4)/(delta_len*del4_len);
        d5=dot(delta(ih,:),del5)/(delta_len*del5_len);
        d6=dot(delta(ih,:),del6)/(delta_len*del6_len);
        if(abs(d4)-1>tol && abs(d5)-1>tol && abs(d6)-1>tol)
            error('Error in build_H.m: hopping cannot be assigned!')
        end
        
        if(sym_sign(ih) < 0)
            % Only hopping of type 4
            if(abs(d4-1)<tol  || abs(d4+1)<tol)
                hopping(ih) = type4(ir,jr(ih),layer);
		drr0 = norm(delta(ih,:) - sign(d4)*del4);
            elseif(abs(d6-1)<tol  || abs(d6+1)<tol)
                hopping(ih) = -type4(ir,jr(ih),layer);
		drr0 = norm(delta(ih,:) - sign(d6)*del6);
            else
                hopping(ih) = 0;
            end
        elseif(sym_sign(ih) > 0)
            % type 4 and type 5
            if(abs(d4-1)<tol || abs(d4+1)<tol)
                hopping(ih) = type4(ir,jr(ih),layer);
		drr0 = norm(delta(ih,:) - sign(d4)*del4);
	    elseif(abs(d6-1)<tol || abs(d6+1)<tol)
                hopping(ih) = type4(ir,jr(ih),layer);
		drr0 = norm(delta(ih,:) - sign(d6)*del6);
            elseif(abs(d5-1)<tol || abs(d5+1)<tol)
                hopping(ih) = type5(ir,jr(ih),layer);
		drr0 = norm(delta(ih,:) - sign(d5)*del5);
            else
                hopping(ih) = 0;
            end
        end
	if(drr0>0.2)
		drr0
		neigh_shell
		error('Something is wrong with drr0')
	end
        hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*drr0/del4_len);
    end
elseif(neigh_shell == 2)
    del1 = a1;
    del2 = a1 + a2;
    del3 = a2;
    del1_len = norm(del1);
    del2_len = norm(del2);
    del3_len = norm(del3);
    for ih = 1 : length(hopping)
	delta_len = norm(delta(ih,:));
        d1=dot(delta(ih,:),del1)/(delta_len*del1_len);
        d2=dot(delta(ih,:),del2)/(delta_len*del2_len);
        d3=dot(delta(ih,:),del3)/(delta_len*del3_len);
        if(abs(d1)-1>tol && abs(d2)-1>tol && abs(d3)-1>tol)
            error('Error in build_H.m: hopping cannot be assigned!')
        end
        if(ir == jr(ih))
            if(abs(d1-1)<tol || abs(d1+1)<tol)
                hopping(ih) = type1(ir,jr(ih),layer);
		drr0 = norm(delta(ih,:) - sign(d1)*del1);
            elseif(abs(d2-1)<tol || abs(d2+1)<tol)
                hopping(ih) = type2(ir,jr(ih),layer);
		drr0 = norm(delta(ih,:) - sign(d2)*del2);
            elseif(abs(d3-1)<tol || abs(d3+1)<tol)
                hopping(ih) = type2(ir,jr(ih),layer);
		drr0 = norm(delta(ih,:) - sign(d3)*del3);
            else
                hopping(ih) = 0;
            end
        elseif(jr(ih) > ir)
            if(sym_sign(ih) < 0)
                %Hopping of type 1 2 and 3
                if(abs(d1+1)<tol)
                    hopping(ih)  = type1(ir,jr(ih),layer);
		    drr0 = norm(delta(ih,:) - sign(d1)*del1);
                elseif(abs(d1-1)<tol)
                    hopping(ih)  = -type1(ir,jr(ih),layer);
		    drr0 = norm(delta(ih,:) - sign(d1)*del1);
                elseif(abs(d2+1)<tol)
                    hopping(ih)  = type2(ir,jr(ih),layer);
		    drr0 = norm(delta(ih,:) - sign(d2)*del2);
                elseif(abs(d3+1)<tol)
                    hopping(ih)  = -type2(ir,jr(ih),layer);
		    drr0 = norm(delta(ih,:) - sign(d3)*del3);
                elseif(abs(d2-1)<tol)
                    hopping(ih)  = -type3(ir,jr(ih),layer);
		    drr0 = norm(delta(ih,:) - sign(d2)*del2);
                elseif(abs(d3-1)<tol)
                    hopping(ih)  = type3(ir,jr(ih),layer);
		    drr0 = norm(delta(ih,:) - sign(d3)*del3);
                else
                    hopping(ih)  = 0;
                end
            elseif(sym_sign(ih) > 0)
                if(abs(d1-1)<tol || abs(d1+1)<tol)
                    hopping(ih) = type1(ir,jr(ih),layer);
		    drr0 = norm(delta(ih,:) - sign(d1)*del1);
                elseif(abs(d2+1)<tol)
                    hopping(ih) = type2(ir,jr(ih),layer);
		    drr0 = norm(delta(ih,:) - sign(d2)*del2);
                elseif(abs(d3+1)<tol)
                    hopping(ih) = type2(ir,jr(ih),layer);
		    drr0 = norm(delta(ih,:) - sign(d3)*del3);
                elseif(abs(d2-1)<tol)
                    hopping(ih) = type3(ir,jr(ih),layer);
		    drr0 = norm(delta(ih,:) - sign(d2)*del2);
                elseif(abs(d3-1)<tol)
                    hopping(ih) = type3(ir,jr(ih),layer);
		    drr0 = norm(delta(ih,:) - sign(d3)*del3);
                else
                    hopping(ih) = 0;
                end
            end
        end
	if(drr0>0.2)
		drr0
		neigh_shell
		error('Something is wrong with drr0')
	end
        hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*drr0/del1_len);
    end
elseif(neigh_shell == 3)
    del7=-2*(a1+2*a2)/3;
    del8=2*(2*a1+a2)/3;
    del9=2*(a2-a1)/3;
    del7_len=norm(del7);
    del8_len=norm(del8);
    del9_len=norm(del9);
    for ih = 1 : length(hopping)
	delta_len = norm(delta(ih,:));
        d7=dot(delta(ih,:),del7)/(delta_len*del7_len);
        d8=dot(delta(ih,:),del8)/(delta_len*del8_len);
        d9=dot(delta(ih,:),del9)/(delta_len*del9_len);
        if(abs(d7)-1>tol && abs(d8)-1>tol && abs(d9)-1>tol)
            error('Error in build_H.m: hopping cannot be assigned!')
        end
        if(sym_sign(ih) < 0)
            % Only hopping of type 8
            if(abs(d8-1)<tol  || abs(d8+1)<tol)
                hopping(ih) = type8(ir,jr(ih),layer);
		drr0 = norm(delta(ih,:) - sign(d8)*del8);
            elseif(abs(d9-1)<tol  || abs(d9+1)<tol)
                hopping(ih) = -type8(ir,jr(ih),layer);
		drr0 = norm(delta(ih,:) - sign(d9)*del9);
            else
                hopping(ih) = 0;
            end
        elseif(sym_sign(ih) > 0)
            % type 8 and type 7
            if(abs(d8-1)<tol  || abs(d8+1)<tol)
                hopping(ih) = type8(ir,jr(ih),layer);
		drr0 = norm(delta(ih,:) - sign(d8)*del8);
	    elseif(abs(d9-1)<tol  || abs(d9+1)<tol)
                hopping(ih) = type8(ir,jr(ih),layer);
		drr0 = norm(delta(ih,:) - sign(d9)*del9);
            elseif(abs(d7-1)<tol || abs(d7+1)<tol)
                hopping(ih) = type7(ir,jr(ih),layer);
		drr0 = norm(delta(ih,:) - sign(d7)*del7);
            else
                hopping(ih) = 0;
            end
        end
        hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*drr0/del7_len);
    end
elseif(neigh_shell == 4)
    tol2 = tol/5;
    del10=-(5*a1+ a2)/3; 
    del11= (a1 + 5*a2)/3; 
    del12= (4*a1 + 5*a2)/3;
    del13= (4*a1 - a2)/3;
    del14= (a1 - 4*a2)/3;
    del15=-(5*a1 + 4*a2)/3;
    del10_len= norm(del10);
    del11_len= norm(del11);
    del12_len= norm(del12);
    del13_len= norm(del13);
    del14_len= norm(del14);
    del15_len= norm(del15);
    for ih = 1 : length(hopping)
	delta_len = norm(delta(ih,:));
        d10=dot(delta(ih,:),del10)/(delta_len*del10_len);
        d11=dot(delta(ih,:),del11)/(delta_len*del11_len);
        d12=dot(delta(ih,:),del12)/(delta_len*del12_len);
        d13=dot(delta(ih,:),del13)/(delta_len*del13_len);
        d14=dot(delta(ih,:),del14)/(delta_len*del14_len);
        d15=dot(delta(ih,:),del15)/(delta_len*del15_len);
        if(abs(d10)-1>tol2 && abs(d11)-1>tol2 && abs(d12)-1>tol2 && ...
		abs(d13)-1>tol2 && abs(d14)-1>tol2 && abs(d15)-1>tol2)
            error('Error in build_H.m: hopping cannot be assigned!')
        end

        if(sym_sign(ih) > 0)
            if(abs(d14-1)<tol2 || abs(d14+1)<tol2)
               hopping(ih) = type14(ir,jr(ih),layer);
	       drr0 = norm(delta(ih,:) - sign(d14)*del14);
	    elseif(abs(d15-1)<tol2 || abs(d15+1)<tol2)
               hopping(ih) = type14(ir,jr(ih),layer);
	       drr0 = norm(delta(ih,:) - sign(d15)*del15);
            elseif(abs(d10-1)<tol2 || abs(d10+1)<tol2)
               hopping(ih) = type10(ir,jr(ih),layer);
	       drr0 = norm(delta(ih,:) - sign(d10)*del10);
	    elseif(abs(d13-1)<tol2 || abs(d13+1)<tol2)
               hopping(ih) = type10(ir,jr(ih),layer);
	       drr0 = norm(delta(ih,:) - sign(d13)*del13);
            elseif(abs(d11-1)<tol2 || abs(d11+1)<tol2)
               hopping(ih) = type12(ir,jr(ih),layer);
	       drr0 = norm(delta(ih,:) - sign(d11)*del11);
	    elseif(abs(d12-1)<tol2 || abs(d12+1)<tol2)
               hopping(ih) = type12(ir,jr(ih),layer);
	       drr0 = norm(delta(ih,:) - sign(d12)*del12);
            else
               hopping(ih) = 0.0;
            end
        elseif(sym_sign(ih) < 0)
            if(abs(d14-1)<tol2 || abs(d14+1)<tol2)
               hopping(ih) = type14(ir,jr(ih),layer);
	       drr0 = norm(delta(ih,:) - sign(d14)*del14);
            elseif(abs(d15-1)<tol2 || abs(d15+1)<tol2)
               hopping(ih) = -type14(ir,jr(ih),layer);
	       drr0 = norm(delta(ih,:) - sign(d15)*del15);
            elseif(abs(d10-1)<tol2 || abs(d10+1)<tol2)
               hopping(ih) = type10(ir,jr(ih),layer);
	       drr0 = norm(delta(ih,:) - sign(d10)*del10);
            elseif(abs(d13-1)<tol2 || abs(d13+1)<tol2)
               hopping(ih) = -type10(ir,jr(ih),layer);
	       drr0 = norm(delta(ih,:) - sign(d13)*del13);
            elseif(abs(d11-1)<tol2 || abs(d11+1)<tol2)
               hopping(ih) = -type12(ir,jr(ih),layer);
	       drr0 = norm(delta(ih,:) - sign(d11)*del11);
            elseif(abs(d12-1)<tol2 || abs(d12+1)<tol2)
               hopping(ih) = type12(ir,jr(ih),layer);
	       drr0 = norm(delta(ih,:) - sign(d12)*del12);
            else
               hopping(ih) = 0.0;
            end
        end
	if(drr0>0.2)
		drr0
		neigh_shell
		error('Something is wrong with drr0')
	end
        hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*drr0/del10_len);
    end
elseif(neigh_shell == 5)
    del16=2*a1 + a2;
    del17=-a1 + a2;
    del18= (2*a2 + a1);
    del16_len=norm(del16);
    del17_len=norm(del17);
    del18_len=norm(del18);
    for ih = 1 : length(hopping)
	delta_len = norm(delta(ih,:));
        d16=dot(delta(ih,:),del16)/(delta_len*del16_len);
        d17=dot(delta(ih,:),del17)/(delta_len*del17_len);
        d18=dot(delta(ih,:),del18)/(delta_len*del18_len);
        if(abs(d16)-1>tol && abs(d17)-1>tol && abs(d18)-1>tol)
            error('Error in build_H.m: hopping cannot be assigned!')
        end
        if(ir == jr(ih))
	   if(abs(d18-1)<tol || abs(d18+1)<tol)
	      hopping(ih) = type18(ir,ir,layer);
	      drr0 = norm(delta(ih,:) - sign(d18)*del18);
           elseif(abs(d16-1)<tol || abs(d16+1)<tol)
	      hopping(ih) = type16(ir,ir,layer);
	      drr0 = norm(delta(ih,:) - sign(d16)*del16);
           elseif(abs(d17-1)<tol || abs(d17+1)<tol)
	      hopping(ih) = type16(ir,ir,layer);
	      drr0 = norm(delta(ih,:) - sign(d17)*del17);
           else
	      hopping(ih) = 0.0;
           end   
        elseif(jr(ih) > ir)
	   if(sym_sign(ih) > 0)
	      if(abs(d16-1)<tol)
	         hopping(ih) = type16(ir,jr(ih),layer); %type16(ir,jr(ih),layer);
	         drr0 = norm(delta(ih,:) - sign(d16)*del16);
	      elseif(abs(d17-1)<tol)
	         hopping(ih) = type16(ir,jr(ih),layer); %type16(ir,jr(ih),layer);
	         drr0 = norm(delta(ih,:) - sign(d17)*del17);
	      elseif(abs(d16+1)<tol)
	         hopping(ih) = typem16(ir,jr(ih),layer); %typem16(ir,jr(ih),layer);
	         drr0 = norm(delta(ih,:) - sign(d16)*del16);
	      elseif(abs(d17+1)<tol)
	         hopping(ih) = typem16(ir,jr(ih),layer); %typem16(ir,jr(ih),layer);
	         drr0 = norm(delta(ih,:) - sign(d17)*del17);
	      elseif(abs(d18-1)<tol)
	         hopping(ih) = type18(ir,jr(ih),layer);
	         drr0 = norm(delta(ih,:) - sign(d18)*del18);
	      elseif(abs(d18+1)<tol)
	         hopping(ih) = typem18(ir,jr(ih),layer); %typem18(ir,jr(ih),layer);
	         drr0 = norm(delta(ih,:) - sign(d18)*del18);
	      else
	         hopping(ih) = 0.0;
              end   
           elseif(sym_sign(ih) < 0)
	      if(abs(d16-1)<tol) 
	         hopping(ih) = type16(ir,jr(ih),layer);
	         drr0 = norm(delta(ih,:) - sign(d16)*del16);
	      elseif(abs(d17-1)<tol)
	         hopping(ih) = -type16(ir,jr(ih),layer);
	         drr0 = norm(delta(ih,:) - sign(d17)*del17);
	      elseif(abs(d16+1)<tol)
	         hopping(ih) = typem16(ir,jr(ih),layer);
	         drr0 = norm(delta(ih,:) - sign(d16)*del16);
              elseif(abs(d17+1)<tol)
	         hopping(ih) = -typem16(ir,jr(ih),layer);
	         drr0 = norm(delta(ih,:) - sign(d17)*del17);
	      else
	         hopping(ih) = 0.0;
              end
           end
        end
	if(drr0>0.2)
		drr0
		neigh_shell
		error('Something is wrong with drr0')
	end
        hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*drr0/del16_len);
     end
elseif(neigh_shell == 6)
    del19= -2*a2; 
    del20= -2*a1; 
    del21= -2*(a1 + a2);
    del19_len=norm(del19);
    del20_len=norm(del20);
    del21_len=norm(del21);
    for ih = 1 : length(hopping)
	delta_len = norm(delta(ih,:));
        d19=dot(delta(ih,:),del19)/(delta_len*del19_len);
        d20=dot(delta(ih,:),del20)/(delta_len*del20_len);
        d21=dot(delta(ih,:),del21)/(delta_len*del21_len);
        if(abs(d19)-1>tol && abs(d20)-1>tol && abs(d21)-1>tol)
            error('Error in build_H.m: hopping cannot be assigned!')
        end
        if(ir == jr(ih))
           if(abs(d20-1)<tol || abs(d20+1)<tol)
	      hopping(ih) = type20(ir,ir,layer);
	      drr0 = norm(delta(ih,:) - sign(d20)*del20);
           elseif(abs(d21-1)<tol || abs(d21+1)<tol)
	      hopping(ih) = type19(ir,ir,layer);
	      drr0 = norm(delta(ih,:) - sign(d21)*del21);
           elseif(abs(d19-1)<tol || abs(d19+1)<tol)
	      hopping(ih) = type19(ir,ir,layer);
	      drr0 = norm(delta(ih,:) - sign(d19)*del19);
           else
	      hopping(ih) = 0.0;
           end
	elseif(jr(ih) > ir)
 	   if(sym_sign(ih) > 0)
              if(abs(d20-1)<tol || abs(d20+1)<tol)
	         hopping(ih) = type20(ir,jr(ih),layer);
	         drr0 = norm(delta(ih,:) - sign(d20)*del20);
	      elseif(abs(d19-1)<tol)
	         hopping(ih) = type19(ir,jr(ih),layer);
	         drr0 = norm(delta(ih,:) - sign(d19)*del19);
	      elseif(abs(d21-1)<tol)
	         hopping(ih) = type19(ir,jr(ih),layer);
	         drr0 = norm(delta(ih,:) - sign(d21)*del21);
	      elseif(abs(d19+1)<tol)
	         hopping(ih) = typem21(ir,jr(ih),layer);
	         drr0 = norm(delta(ih,:) - sign(d19)*del19);
	      elseif(abs(d21+1)<tol)
	         hopping(ih) = typem21(ir,jr(ih),layer);
	         drr0 = norm(delta(ih,:) - sign(d21)*del21);
	      else
		 hopping(ih) = 0.0;
	      end	 
	   elseif(sym_sign(ih) < 0)
              if(abs(d20+1)<tol) 
	         hopping(ih) = -type20(ir,jr(ih),layer);
	         drr0 = norm(delta(ih,:) - sign(d20)*del20);
	      elseif(abs(d20-1)<tol) 
	         hopping(ih) = type20(ir,jr(ih),layer);
	         drr0 = norm(delta(ih,:) - sign(d20)*del20);
	      elseif(abs(d19-1)<tol)
	         hopping(ih) = type19(ir,jr(ih),layer);
	         drr0 = norm(delta(ih,:) - sign(d19)*del19);
	      elseif(abs(d21-1)<tol)
	         hopping(ih) = -type19(ir,jr(ih),layer);
	         drr0 = norm(delta(ih,:) - sign(d21)*del21);
	      elseif(abs(d19+1)<tol) 
	         hopping(ih) = -typem21(ir,jr(ih),layer);
	         drr0 = norm(delta(ih,:) - sign(d19)*del19);
	      elseif(abs(d21+1)<tol)
	         hopping(ih) = typem21(ir,jr(ih),layer);
	         drr0 = norm(delta(ih,:) - sign(d21)*del21);
	      else
	         hopping(ih) = 0.0;
	      end
           end
       end
	if(drr0>0.2)
		drr0
		neigh_shell
		error('Something is wrong with drr0')
	end
       hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*drr0/del19_len);
    end
end     
end

