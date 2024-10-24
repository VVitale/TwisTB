function hopping = assign_hoppings(a1,a2,delta,loc_theta1,loc_theta2,neigh_shell,layer,...
    sym_sign,ir,jr,type1,type2,type3,type4,type5,type7,type8,type10,type12,...
    type14,type16,typem16,type18,typem18,type19,type20,typem21)
hopping = zeros(size(delta,1),1);
tol = 0.3;
c2t1 = cos(2*loc_theta1);
s2t1 = sin(2*loc_theta1);
c2t2 = cos(2*loc_theta2);
s2t2 = sin(2*loc_theta2);
if (neigh_shell == 1)
    del4=-(2*a1+a2)/3;
    del5=(a1+2*a2)/3;
    del6=(a1-a2)/3;
    for ih = 1 : length(hopping)
        d4=dot(delta(ih,:),del4)/(norm(delta(ih,:))*norm(del4));
        d5=dot(delta(ih,:),del5)/(norm(delta(ih,:))*norm(del5));
        d6=dot(delta(ih,:),del6)/(norm(delta(ih,:))*norm(del6));
        if(abs(d4)-1>tol && abs(d5)-1>tol && abs(d6)-1>tol)
            error('Error in build_H.m: hopping cannot be assigned!')
        end
        
        if(sym_sign(ih) < 0)
            % Only hopping of type 4
            if(abs(d4-1)<tol  || abs(d4+1)<tol)
                hopping(ih) = type4(ir,jr(ih),layer);
                if(ir == 9 && jr(ih) == 7)
                    hopping(ih) = type4(9,7,layer)*c2t2(ih) - type4(9,8)*s2t2(ih);
                elseif(ir == 7 && jr(ih) == 9)
                    hopping(ih) = type4(7,9,layer)*c2t1 - type4(8,9)*s2t1;
                elseif(ir == 11 && jr(ih) == 7)
                    hopping(ih) = type4(11,7,layer)*c2t2(ih) - type4(11,8)*s2t2(ih);
                elseif(ir == 7 && jr(ih) == 11)
                    hopping(ih) = type4(7,11,layer)*c2t1 - type4(8,11)*s2t1;
                elseif(ir == 10 && jr(ih) == 8)
                    hopping(ih) = type4(10,8,layer)*c2t2(ih) + type4(10,7)*s2t2(ih);
                elseif(ir == 8  && jr(ih) == 10)
                    hopping(ih) = type4(8,10,layer)*c2t1 + type4(7,10)*s2t1;
                end
                hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*(norm(delta(ih,:)) - norm(del4))/norm(del4));
            elseif(abs(d6-1)<tol  || abs(d6+1)<tol)
                hopping(ih) = -type4(ir,jr(ih),layer);
                if(ir == 9 && jr(ih) == 7)
                    hopping(ih) = -type4(9,7,layer)*c2t2(ih) - type4(9,8)*s2t2(ih);
                elseif(ir == 7 && jr(ih) == 9)
                    hopping(ih) = -type4(7,9,layer)*c2t1 - type4(8,9)*s2t1;
                elseif(ir == 11 && jr(ih) == 7)
                    hopping(ih) = -type4(11,7,layer)*c2t2(ih) - type4(11,8)*s2t2(ih);
                elseif(ir == 7 && jr(ih) == 11)
                    hopping(ih) = -type4(7,11,layer)*c2t1 - type4(8,11)*s2t1;
                elseif(ir == 10 && jr(ih) == 8)
                    hopping(ih) = -type4(10,8,layer)*c2t2(ih) + type4(10,7)*s2t2(ih);
                elseif(ir == 8 && jr(ih) == 10)
                    hopping(ih) = -type4(8,10,layer)*c2t1 + type4(7,10)*s2t1;
                end
                hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*(norm(delta(ih,:)) - norm(del6))/norm(del6));
            elseif(abs(d5-1) < tol || abs(d5+1) < tol)
                if(ir == 9 && jr(ih) ==7)
                    hopping(ih) = -type5(9,8,layer)*s2t2(ih);
                elseif(ir == 7 && jr(ih) == 9)
                    hopping(ih) = -type5(8,9,layer)*s2t1;
                elseif(ir == 11 && jr(ih) ==7)
                    hopping(ih) = -type5(11,8,layer)*s2t2(ih);
                elseif(ir == 7  && jr(ih) ==11)
                    hopping(ih) = -type5(8,11,layer)*s2t1;
                elseif(ir == 10 && jr(ih) ==8)
                    hopping(ih) = type5(10,7,layer)*s2t2(ih);
                elseif(ir == 8  && jr(ih) ==10)
                    hopping(ih) = type5(7,10,layer)*s2t1;
                end
                hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*(norm(delta(ih,:)) - norm(del5))/norm(del5));
            else
                hopping(ih) = 0;
            end
        elseif(sym_sign(ih) > 0)
            % type 4 and type 5
            if(abs(d4-1)<tol || abs(d4+1)<tol)
                hopping(ih) = type4(ir,jr(ih),layer);
                if(ir == 10 && jr(ih) == 7)
                    hopping(ih) = type4(10,7,layer)*c2t2(ih) - type4(10,8)*s2t2(ih);
                elseif(ir == 7 && jr(ih) == 10)
                    hopping(ih) = type4(7,10,layer)*c2t1 - type4(8,10)*s2t1;
                elseif(ir == 9 && jr(ih) == 8)
                    hopping(ih) = type4(9,8,layer)*c2t2(ih) + type4(9,7)*s2t2(ih);
                elseif(ir == 8 && jr(ih) == 9)
                    hopping(ih) = type4(8,9,layer)*c2t1 + type4(7,9)*s2t1;
                elseif(ir == 11 && jr(ih) == 8)
                    hopping(ih) = type4(11,8,layer)*c2t2(ih) + type4(11,7)*s2t2(ih);
                elseif(ir == 8 && jr(ih) == 11)
                    hopping(ih) = type4(8,11,layer)*c2t1 + type4(7,11)*s2t1;
                end
                hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*(norm(delta(ih,:)) - norm( del4))/norm(del4));
            elseif(abs(d6-1)<tol || abs(d6+1)<tol)
                hopping(ih) = type4(ir,jr(ih),layer);
                if(ir == 10 && jr(ih) == 7)
                    hopping(ih) = type4(10,7,layer)*c2t2(ih) + type4(10,8)*s2t2(ih);
                elseif(ir == 7 && jr(ih) == 10)
                    hopping(ih) = type4(7,10,layer)*c2t1 + type4(8,10)*s2t1;
                elseif(ir == 9 && jr(ih) == 8)
                    hopping(ih) = type4(9,8,layer)*c2t2(ih) - type4(9,7)*s2t2(ih);
                elseif(ir == 8 && jr(ih) == 9)
                    hopping(ih) = type4(8,9,layer)*c2t1 - type4(7,9)*s2t1;
                elseif(ir == 11 && jr(ih) == 8)
                    hopping(ih) = type4(11,8,layer)*c2t2(ih) - type4(11,7)*s2t2(ih);
                elseif(ir == 8 && jr(ih) == 11)
                    hopping(ih) = type4(8,11,layer)*c2t1 - type4(7,11)*s2t1;
                end
                hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*(norm(delta(ih,:)) - norm( del6))/norm(del6));
            elseif(abs(d5-1)<tol || abs(d5+1)<tol)
                hopping(ih) = type5(ir,jr(ih),layer);
                if(ir == 10 && jr(ih) == 7)
                    hopping(ih) = type5(10,7,layer)*c2t2(ih);
                elseif(ir == 7 && jr(ih) == 10)
                    hopping(ih) = type5(7,10,layer)*c2t1;
                elseif(ir == 9 && jr(ih) == 8)
                    hopping(ih) = type5(9,8,layer)*c2t2(ih);
                elseif(ir == 8 && jr(ih) == 9)
                    hopping(ih) = type5(8,9,layer)*c2t1;
                elseif(ir == 11 && jr(ih) == 8)
                    hopping(ih) = type5(11,8,layer)*c2t2(ih);
                elseif(ir == 8 && jr(ih) == 11)
                    hopping(ih) = type5(8,11,layer)*c2t1;
                end
                hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*(norm(delta(ih,:)) - norm(del5))/norm(del5));
            else
                hopping(ih) = 0;
            end
        end
    end
elseif(neigh_shell == 2)
    del1 = a1;
    del2 = a1 + a2;
    del3 = a2;
    for ih = 1 : length(hopping)
        d1=dot(delta(ih,:),del1)/(norm(delta(ih,:))*norm(del1));
        d2=dot(delta(ih,:),del2)/(norm(delta(ih,:))*norm(del2));
        d3=dot(delta(ih,:),del3)/(norm(delta(ih,:))*norm(del3));
        if(abs(d1)-1>tol && abs(d2)-1>tol && abs(d3)-1>tol)
            error('Error in build_H.m: hopping cannot be assigned!')
        end
        if(ir == jr(ih))
            if(abs(d1-1)<tol)
                hopping(ih) = type1(ir,jr(ih),layer);
                if(ir == 7)
                    hopping(ih) = type1(7,7,layer)*c2t1*c2t2(ih) + type1(8,8,layer)*s2t1*s2t2(ih)  - type1(7,8,layer)*c2t1*s2t2(ih) + type1(8,7,layer)*s2t1*c2t2(ih);
                elseif(ir == 8)
                    hopping(ih) = type1(8,8,layer)*c2t1*c2t2(ih) + type1(7,7,layer)*s2t1*s2t2(ih)  - type1(8,7,layer)*c2t1*s2t2(ih) + type1(7,8,layer)*s2t1*c2t2(ih);
                end
                hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*(norm(delta(ih,:)) - norm(del1))/norm(del1));
            elseif(abs(d1+1)<tol)
                hopping(ih) = type1(ir,jr(ih),layer);
                if(ir == 7)
                    hopping(ih) = type1(7,7,layer)*c2t1*c2t2(ih) + type1(8,8,layer)*s2t1*s2t2(ih)  + type1(7,8,layer)*c2t1*s2t2(ih) - type1(8,7,layer)*s2t1*c2t2(ih);
                elseif(ir == 8)
                    hopping(ih) = type1(8,8,layer)*c2t1*c2t2(ih) + type1(7,7,layer)*s2t1*s2t2(ih)  + type1(8,7,layer)*c2t1*s2t2(ih) - type1(7,8,layer)*s2t1*c2t2(ih);
                end
                hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*(norm(delta(ih,:)) - norm(del1))/norm(del1));
            elseif(abs(d2-1)<tol)
                hopping(ih) = type2(ir,jr(ih),layer);
                if(ir == 7)
                    hopping(ih) = type2(7,7,layer)*c2t1*c2t2(ih) + type2(8,8,layer)*s2t1*s2t2(ih) + type3(7,8,layer)*c2t1*s2t2(ih) - type2(8,7,layer)*s2t1*c2t2(ih);
                elseif(ir == 8)
                    hopping(ih) = type2(8,8,layer)*c2t1*c2t2(ih) + type2(7,7,layer)*s2t1*s2t2(ih) + type2(8,7,layer)*c2t1*s2t2(ih) - type3(7,8,layer)*s2t1*c2t2(ih);
                end
                hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*(norm(delta(ih,:)) - norm(del2))/norm(del2));
            elseif(abs(d2+1)<tol)
                hopping(ih) = type2(ir,jr(ih),layer);
                if(ir == 7)
                    hopping(ih) = type2(7,7,layer)*c2t1*c2t2(ih) + type2(8,8,layer)*s2t1*s2t2(ih) - type2(7,8,layer)*c2t1*s2t2(ih) + type3(8,7,layer)*s2t1*c2t2(ih);
                elseif(ir == 8)
                    hopping(ih) = type2(8,8,layer)*c2t1*c2t2(ih) + type2(7,7,layer)*s2t1*s2t2(ih) - type3(8,7,layer)*c2t1*s2t2(ih) + type2(7,8,layer)*s2t1*c2t2(ih);
                end
                hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*(norm(delta(ih,:)) - norm(del2))/norm(del2));
            elseif(abs(d3-1)<tol)
                hopping(ih) = type2(ir,jr(ih),layer);
                if(ir == 7)
                    hopping(ih) = type2(7,7,layer)*c2t1*c2t2(ih) + type2(8,8,layer)*s2t1*s2t2(ih) - type3(7,8,layer)*c2t1*s2t2(ih) + type2(8,7,layer)*s2t1*c2t2(ih);
                elseif(ir == 8)
                    hopping(ih) = type2(8,8,layer)*c2t1*c2t2(ih) + type2(7,7,layer)*s2t1*s2t2(ih) - type2(8,7,layer)*c2t1*s2t2(ih) + type3(7,8,layer)*s2t1*c2t2(ih);
                end
                hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*(norm(delta(ih,:)) - norm(del3))/norm(del3));
            elseif(abs(d3+1)<tol)
                hopping(ih) = type2(ir,jr(ih),layer);
                if(ir == 7)
                    hopping(ih) = type2(7,7,layer)*c2t1*c2t2(ih) + type2(8,8,layer)*s2t1*s2t2(ih) + type2(7,8,layer)*c2t1*s2t2(ih) - type3(8,7,layer)*s2t1*c2t2(ih);
                elseif(ir == 8)
                    hopping(ih) = type2(8,8,layer)*c2t1*c2t2(ih) + type2(7,7,layer)*s2t1*s2t2(ih) + type3(8,7,layer)*c2t1*s2t2(ih) - type2(7,8,layer)*s2t1*c2t2(ih);
                end
                hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*(norm(delta(ih,:)) - norm(del3))/norm(del3));
            else
                hopping(ih) = 0;
            end
        elseif(jr(ih) > ir)
            if(sym_sign(ih) < 0)
                %Hopping of type 1 2 and 3
                if(abs(d1-1)<tol)
                    hopping(ih)  = -type1(ir,jr(ih),layer);
                    if(jr(ih)==8 && ir == 7)
                        hopping(ih) = -type1(7,8,layer)*c2t1*c2t2(ih) + type1(8,7,layer)*s2t1*s2t2(ih) + type1(7,7,layer)*c2t1*s2t2(ih) - type1(8,8,layer)*s2t1*c2t2(ih);
                    elseif(jr(ih)==7 && ir == 6)
                        hopping(ih) = -type1(6,7,layer)*c2t2(ih) - type1(6,8,layer)*s2t2(ih);
                    end
                    hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*(norm(delta(ih,:)) - norm(del1))/norm(del1));
                elseif(abs(d1+1)<tol)
                    hopping(ih)  = type1(ir,jr(ih),layer);
                    if(jr(ih)==8 && ir == 7)
                        hopping(ih) =  type1(7,8,layer)*c2t1*c2t2(ih) - type1(8,7,layer)*s2t1*s2t2(ih) + type1(7,7,layer)*c2t1*s2t2(ih) - type1(8,8,layer)*s2t1*c2t2(ih);
                    elseif(jr(ih)==7 && ir == 6)
                        hopping(ih) = type1(6,7,layer)*c2t2(ih) - type1(6,8,layer)*s2t2(ih);
                    end
                    hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*(norm(delta(ih,:)) - norm(del1))/norm(del1));
                elseif(abs(d2+1)<tol)
                    hopping(ih)  = type2(ir,jr(ih),layer);
                    if(jr(ih)==8 && ir == 7)
                        hopping(ih) = type2(7,8,layer)*c2t1*c2t2(ih) - type2(8,7,layer)*s2t1*s2t2(ih) + type2(7,7,layer)*c2t1*s2t2(ih) - type2(8,8,layer)*s2t1*c2t2(ih);
                    elseif(jr(ih)==7 && ir == 6)
                        hopping(ih) = type2(6,7,layer)*c2t2(ih) - type2(6,8,layer)*s2t2(ih);
                    end
                    hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*(norm(delta(ih,:)) - norm(del2))/norm(del2));
                elseif(abs(d3+1)<tol)
                    hopping(ih)  = -type2(ir,jr(ih),layer);
                    if(jr(ih)==8 && ir == 7)
                        hopping(ih) = -type2(7,8,layer)*c2t1*c2t2(ih) + type2(8,7,layer)*s2t1*s2t2(ih) + type2(7,7,layer)*c2t1*s2t2(ih) - type2(8,8,layer)*s2t1*c2t2(ih);
                    elseif(jr(ih)==7 && ir == 6)
                        hopping(ih) = -type2(6,7,layer)*c2t2(ih) - type2(6,8,layer)*s2t2(ih);
                    end
                    hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*(norm(delta(ih,:)) - norm(del3))/norm(del3));
                elseif(abs(d2-1)<tol)
                    hopping(ih)  = -type3(ir,jr(ih),layer);
                    if(jr(ih)==8 && ir == 7)
                        hopping(ih) = -type3(7,8,layer)*c2t1*c2t2(ih) + type3(8,7,layer)*s2t1*s2t2(ih) + type2(7,7,layer)*c2t1*s2t2(ih) - type2(8,8,layer)*s2t1*c2t2(ih);
                    elseif(jr(ih)==7 && ir == 6)
                        hopping(ih) = -type3(6,7,layer)*c2t2(ih) - type3(6,8,layer)*s2t2(ih);
                    end
                    hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*(norm(delta(ih,:)) - norm(del2))/norm(del2));
                elseif(abs(d3-1)<tol)
                    hopping(ih)  = type3(ir,jr(ih),layer);
                    if(jr(ih)==8 && ir == 7)
                        hopping(ih) = type3(7,8,layer)*c2t1*c2t2(ih) - type3(8,7,layer)*s2t1*s2t2(ih) + type2(7,7,layer)*c2t1*s2t2(ih) - type2(8,8,layer)*s2t1*c2t2(ih);
                    elseif(jr(ih)==7 && ir == 6)
                        hopping(ih) = type3(6,7,layer)*c2t2(ih) - type3(6,8,layer)*s2t2(ih);
                    end
                    hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*(norm(delta(ih,:)) - norm(del3))/norm(del3));
                else
                    hopping(ih)  = 0;
                end
            elseif(sym_sign(ih) > 0)
                if(abs(d1-1)<tol)
                    hopping(ih) = type1(ir,jr(ih),layer);
                    if(jr(ih)==8 && ir == 6)
                        hopping(ih) = type1(6,8,layer)*c2t2(ih) - type1(6,7,layer)*s2t2(ih);
                    end
                    hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*(norm(delta(ih,:)) - norm(del1))/norm(del1));
                elseif(abs(d1+1)<tol)
                    hopping(ih) = type1(ir,jr(ih),layer);
                    if(jr(ih)==8 && ir == 6)
                        hopping(ih) = type1(6,8,layer)*c2t2(ih) + type1(6,7,layer)*s2t2(ih);
                    end
                    hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*(norm(delta(ih,:)) - norm(del1))/norm(del1));
                elseif(abs(d2+1)<tol)
                    hopping(ih) = type2(ir,jr(ih),layer);
                    if(jr(ih)==8 && ir == 6)
                        hopping(ih) = type2(6,8,layer)*c2t2(ih) + type2(6,7,layer)*s2t2(ih);
                    end
                    hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*(norm(delta(ih,:)) - norm(del2))/norm(del2));
                elseif(abs(d3+1)<tol)
                    hopping(ih) = type2(ir,jr(ih),layer);
                    if(jr(ih)==8 && ir == 6)
                        hopping(ih) = type2(6,8,layer)*c2t2(ih) - type2(6,7,layer)*s2t2(ih);
                    end
                    hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*(norm(delta(ih,:)) - norm(del3))/norm(del3));
                elseif(abs(d2-1)<tol)
                    hopping(ih) = type3(ir,jr(ih),layer);
                    if(jr(ih)==8 && ir == 6)
                        hopping(ih) = type3(6,8,layer)*c2t2(ih) - type3(6,7,layer)*s2t2(ih);
                    end
                    hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*(norm(delta(ih,:)) - norm(del2))/norm(del2));
                elseif(abs(d3-1)<tol)
                    hopping(ih) = type3(ir,jr(ih),layer);
                    if(jr(ih)==8 && ir == 6)
                        hopping(ih) = type3(6,8,layer)*c2t2(ih) + type3(6,7,layer)*s2t2(ih);
                    end
                    hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*(norm(delta(ih,:)) - norm(del3))/norm(del3));
                else
                    hopping(ih) = 0;
                end
            end
        end
    end
elseif(neigh_shell == 3)
    del7=-2*(a1+2*a2)/3;
    del8=2*(2*a1+a2)/3;
    del9=2*(a2-a1)/3;
    for ih = 1 : length(hopping)
        d7=dot(delta(ih,:),del7)/(norm(delta(ih,:))*norm(del7));
        d8=dot(delta(ih,:),del8)/(norm(delta(ih,:))*norm(del8));
        d9=dot(delta(ih,:),del9)/(norm(delta(ih,:))*norm(del9));
        if(abs(d7)-1>tol && abs(d8)-1>tol && abs(d9)-1>tol)
            error('Error in build_H.m: hopping cannot be assigned!')
        end
        if(sym_sign(ih) < 0)
            % Only hopping of type 8
            if(abs(d8-1)<tol  || abs(d8+1)<tol)
                hopping(ih) = type8(ir,jr(ih),layer);
		if(ir == 9 && jr(ih) == 7)
		   hopping(ih) = type8(9,7,layer)*c2t2(ih) - type8(9,8)*s2t2(ih);
	        elseif(ir == 7 && jr(ih) == 9)
		   hopping(ih) = type8(7,9,layer)*c2t1 - type8(8,9)*s2t1; 
	        elseif(ir == 11 && jr(ih) == 7)
		   hopping(ih) = type8(11,7,layer)*c2t2(ih) - type8(11,8)*s2t2(ih); 
	        elseif(ir == 7 && jr(ih) == 11)
		   hopping(ih) = type8(7,11,layer)*c2t1 - type8(8,11)*s2t1;
	        elseif(ir == 10 && jr(ih) == 8)
		   hopping(ih) = type8(10,8,layer)*c2t2(ih) + type8(10,7)*s2t2(ih); 
	        elseif(ir == 8  && jr(ih) == 10)
		   hopping(ih) = type8(8,10,layer)*c2t1 + type8(7,10)*s2t1;
	        end
                hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*(norm(delta(ih,:))- norm(del8))/norm(del8));
            elseif(abs(d9-1)<tol  || abs(d9+1)<tol)
                hopping(ih) = -type8(ir,jr(ih),layer);
		if(ir == 9 && jr(ih) == 7)
		   hopping(ih) = -type8(9,7,layer)*c2t2(ih) - type8(9,8)*s2t2(ih); 
	        elseif(ir == 7 && jr(ih) == 9)
		   hopping(ih) = -type8(7,9,layer)*c2t1 - type8(8,9)*s2t1; 
	        elseif(ir == 11 && jr(ih) == 7)
		   hopping(ih) = -type8(11,7,layer)*c2t2(ih) - type8(11,8)*s2t2(ih); 
	        elseif(ir == 7 && jr(ih) == 11)
		   hopping(ih) = -type8(7,11,layer)*c2t1 - type8(8,11)*s2t1;
	        elseif(ir == 10 && jr(ih) == 8)
		   hopping(ih) = -type8(10,8,layer)*c2t2(ih) + type8(10,7)*s2t2(ih); 
	        elseif(ir == 8 && jr(ih) == 10)
		   hopping(ih) = -type8(8,10,layer)*c2t1 + type8(7,10)*s2t1;
	        end
                hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*(norm(delta(ih,:))- norm(del9))/norm(del9));
	    elseif(abs(d7-1) < tol || abs(d7+1) < tol)
		if(ir == 9 && jr(ih) ==7)
		   hopping(ih) = -type7(9,8,layer)*s2t2(ih);
		elseif(ir == 7 && jr(ih) == 9)
		   hopping(ih) = -type7(8,9,layer)*s2t1;
	        elseif(ir == 11 && jr(ih) ==7)
		   hopping(ih) = -type7(11,8,layer)*s2t2(ih);
	        elseif(ir == 7  && jr(ih) ==11)
		   hopping(ih) = -type7(8,11,layer)*s2t1;
	        elseif(ir == 10 && jr(ih) ==8)
		   hopping(ih) = type7(10,7,layer)*s2t2(ih);
	        elseif(ir == 8  && jr(ih) ==10)
		   hopping(ih) = type7(7,10,layer)*s2t1;
	        end
                hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*(norm(delta(ih,:))- norm(del7))/norm(del7));
            else
                hopping(ih) = 0;
            end
        elseif(sym_sign(ih) > 0)
            % type 8 and type 7
            if(abs(d8-1)<tol || abs(d8+1)<tol)
                hopping(ih) = type8(ir,jr(ih),layer);
		if(ir == 10 && jr(ih) == 7)
		   hopping(ih) = type8(10,7,layer)*c2t2(ih) - type8(10,8)*s2t2(ih); 
	        elseif(ir == 7 && jr(ih) == 10)
		   hopping(ih) = type8(7,10,layer)*c2t1 - type8(8,10)*s2t1; 
	        elseif(ir == 9 && jr(ih) == 8)
		   hopping(ih) = type8(9,8,layer)*c2t2(ih) + type8(9,7)*s2t2(ih); 
	        elseif(ir == 8 && jr(ih) == 9)
		   hopping(ih) = type8(8,9,layer)*c2t1 + type8(7,9)*s2t1; 
	        elseif(ir == 11 && jr(ih) == 8)
		   hopping(ih) = type8(11,8,layer)*c2t2(ih) + type8(11,7)*s2t2(ih); 
	        elseif(ir == 8 && jr(ih) == 11)
		   hopping(ih) = type8(8,11,layer)*c2t1 + type8(7,11)*s2t1;
	        end
                hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*(norm(delta(ih,:))- norm(del8))/norm(del8));
	    elseif(abs(d9-1)<tol || abs(d9+1)<tol)
                hopping(ih) = type8(ir,jr(ih),layer);
		if(ir == 10 && jr(ih) == 7)
		   hopping(ih) = type8(10,7,layer)*c2t2(ih) + type8(10,8)*s2t2(ih); 
	        elseif(ir == 7 && jr(ih) == 10)
		   hopping(ih) = type8(7,10,layer)*c2t1 + type8(8,10)*s2t1; 
	        elseif(ir == 9 && jr(ih) == 8)
		   hopping(ih) = type8(9,8,layer)*c2t2(ih) - type8(9,7)*s2t2(ih); 
	        elseif(ir == 8 && jr(ih) == 9)
		   hopping(ih) = type8(8,9,layer)*c2t1 - type8(7,9)*s2t1; 
	        elseif(ir == 11 && jr(ih) == 8)
		   hopping(ih) = type8(11,8,layer)*c2t2(ih) - type8(11,7)*s2t2(ih); 
	        elseif(ir == 8 && jr(ih) == 11)
		   hopping(ih) = type8(8,11,layer)*c2t1 - type8(7,11)*s2t1;
	        end
                hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*(norm(delta(ih,:))- norm(del9))/norm(del9));
            elseif(abs(d7-1)<tol || abs(d7+1)<tol)
                hopping(ih) = type7(ir,jr(ih),layer);
		if(ir == 10 && jr(ih) == 7)
		   hopping(ih) = type7(10,7,layer)*c2t2(ih); 
	        elseif(ir == 7 && jr(ih) == 10)
		   hopping(ih) = type7(7,10,layer)*c2t1; 
	        elseif(ir == 9 && jr(ih) == 8)
		   hopping(ih) = type7(9,8,layer)*c2t2(ih); 
	        elseif(ir == 8 && jr(ih) == 9)
		   hopping(ih) = type7(8,9,layer)*c2t1; 
	        elseif(ir == 11 && jr(ih) == 8)
		   hopping(ih) = type7(11,8,layer)*c2t2(ih); 
	        elseif(ir == 8 && jr(ih) == 11)
		   hopping(ih) = type7(8,11,layer)*c2t1;
	        end
                hopping(ih) = hopping(ih)*(1.0 - Lambda(ir,jr(ih))*(norm(delta(ih,:))- norm(del7))/norm(del7));
            else
                hopping(ih) = 0;
            end
        end
    end
elseif(neigh_shell == 4)
    tol2 = tol/5;
    del10=-(5*a1+ a2)/3; 
    del11= (a1 + 5*a2)/3; 
    del12= (4*a1 + 5*a2)/3;
    del13= (4*a1 - a2)/3;
    del14= (a1 - 4*a2)/3;
    del15=-(5*a1 + 4*a2)/3;
    for ih = 1 : length(hopping)
        d10=dot(delta(ih,:),del10)/(norm(delta(ih,:))*norm(del10));
        d11=dot(delta(ih,:),del11)/(norm(delta(ih,:))*norm(del11));
        d12=dot(delta(ih,:),del12)/(norm(delta(ih,:))*norm(del12));
        d13=dot(delta(ih,:),del13)/(norm(delta(ih,:))*norm(del13));
        d14=dot(delta(ih,:),del14)/(norm(delta(ih,:))*norm(del14));
        d15=dot(delta(ih,:),del15)/(norm(delta(ih,:))*norm(del15));
        if(abs(d10)-1>tol2 && abs(d11)-1>tol2 && abs(d12)-1>tol2 && ...
		abs(d13)-1>tol2 && abs(d14)-1>tol2 && abs(d15)-1>tol2)
            error('Error in build_H.m: hopping cannot be assigned!')
        end

        if(sym_sign(ih) > 0)
            if(abs(d14-1)<tol2 || abs(d15-1)<tol2 || abs(d14+1)<tol2 || abs(d15+1)<tol2)
               hopping(ih) = type14(ir,jr(ih),layer);
            elseif(abs(d10-1)<tol2 || abs(d13-1)<tol2 || abs(d10+1)<tol2 || abs(d13+1)<tol2)
               hopping(ih) = type10(ir,jr(ih),layer);
            elseif(abs(d11-1)<tol2 || abs(d12-1)<tol2 || abs(d11+1)<tol2 || abs(d12+1)<tol2)
               hopping(ih) = type12(ir,jr(ih),layer);
            else
               hopping(ih) = 0.0;
            end
        elseif(sym_sign(ih) < 0)
            if(abs(d14-1)<tol2 || abs(d14+1)<tol2)
               hopping(ih) = type14(ir,jr(ih),layer);
            elseif(abs(d15-1)<tol2 || abs(d15+1)<tol2)
               hopping(ih) = -type14(ir,jr(ih),layer);
            elseif(abs(d10-1)<tol2 || abs(d10+1)<tol2)
               hopping(ih) = type10(ir,jr(ih),layer);
            elseif(abs(d13-1)<tol2 || abs(d13+1)<tol2)
               hopping(ih) = -type10(ir,jr(ih),layer);
            elseif(abs(d11-1)<tol2 || abs(d11+1)<tol2)
               hopping(ih) = -type12(ir,jr(ih),layer);
            elseif(abs(d12-1)<tol2 || abs(d12+1)<tol2)
               hopping(ih) = type12(ir,jr(ih),layer);
            else
               hopping(ih) = 0.0;
            end
        end
    end
elseif(neigh_shell == 5)
    del16=2*a1 + a2;
    del17=-a1 + a2;
    del18= (2*a2 + a1);
    for ih = 1 : length(hopping)
        d16=dot(delta(ih,:),del16)/(norm(delta(ih,:))*norm(del16));
        d17=dot(delta(ih,:),del17)/(norm(delta(ih,:))*norm(del17));
        d18=dot(delta(ih,:),del18)/(norm(delta(ih,:))*norm(del18));
        if(abs(d16)-1>tol && abs(d17)-1>tol && abs(d18)-1>tol)
            error('Error in build_H.m: hopping cannot be assigned!')
        end
        if(ir == jr(ih))
	   if(abs(d18-1)<tol || abs(d18+1)<tol)
	      hopping(ih) = type18(ir,ir,layer);
           elseif(abs(d16-1)<tol || abs(d16+1)<tol)
	      hopping(ih) = type16(ir,ir,layer);
           elseif(abs(d17-1)<tol || abs(d17+1)<tol)
	      hopping(ih) = type16(ir,ir,layer);
           else
	      hopping(ih) = 0.0;
           end   
        elseif(jr(ih) > ir)
	   if(sym_sign(ih) > 0)
	      if(abs(d16-1)<tol || abs(d17-1)<tol)
	         hopping(ih) = type16(ir,jr(ih),layer); %type16(ir,jr(ih),layer);
	      elseif(abs(d16+1)<tol || abs(d17+1)<tol)
	         hopping(ih) = typem16(ir,jr(ih),layer); %typem16(ir,jr(ih),layer);
	      elseif(abs(d18-1)<tol)
	         hopping(ih) = type18(ir,jr(ih),layer);
	      elseif(abs(d18+1)<tol)
	         hopping(ih) = typem18(ir,jr(ih),layer); %typem18(ir,jr(ih),layer);
	      else
	         hopping(ih) = 0.0;
              end   
           elseif(sym_sign(ih) < 0)
	      if(abs(d16-1)<tol) 
	         hopping(ih) = type16(ir,jr(ih),layer);
	      elseif(abs(d17-1)<tol)
	         hopping(ih) = -type16(ir,jr(ih),layer);
	      elseif(abs(d16+1)<tol)
	         hopping(ih) = typem16(ir,jr(ih),layer);
              elseif(abs(d17+1)<tol)
	         hopping(ih) = -typem16(ir,jr(ih),layer);
	      else
	         hopping(ih) = 0.0;
              end
           end
        end
     end
elseif(neigh_shell == 6)
    del19= -2*a2; 
    del20= -2*a1; 
    del21= -2*(a1 + a2);
    for ih = 1 : length(hopping)
        d19=dot(delta(ih,:),del19)/(norm(delta(ih,:))*norm(del19));
        d20=dot(delta(ih,:),del20)/(norm(delta(ih,:))*norm(del20));
        d21=dot(delta(ih,:),del21)/(norm(delta(ih,:))*norm(del21));
        if(abs(d19)-1>tol && abs(d20)-1>tol && abs(d21)-1>tol)
            error('Error in build_H.m: hopping cannot be assigned!')
        end
        if(ir == jr(ih))
           if(abs(d20-1)<tol || abs(d20+1)<tol)
	      hopping(ih) = type20(ir,ir,layer);
           elseif(abs(d21-1)<tol || abs(d21+1)<tol)
	      hopping(ih) = type19(ir,ir,layer);
           elseif(abs(d19-1)<tol || abs(d19+1)<tol)
	      hopping(ih) = type19(ir,ir,layer);
           else
	      hopping(ih) = 0.0;
           end
	elseif(jr(ih) > ir)
 	   if(sym_sign(ih) > 0)
              if(abs(d20-1)<tol || abs(d20+1)<tol)
	         hopping(ih) = type20(ir,jr(ih),layer);
	      elseif(abs(d19-1)<tol || abs(d21-1)<tol)
	         hopping(ih) = type19(ir,jr(ih),layer);
	      elseif(abs(d19+1)<tol || abs(d21+1)<tol)
	         hopping(ih) = typem21(ir,jr(ih),layer);
	      else
		 hopping(ih) = 0.0;
	      end	 
	   elseif(sym_sign(ih) < 0)
              if(abs(d20+1)<tol) 
	         hopping(ih) = -type20(ir,jr(ih),layer);
	      elseif(abs(d20-1)<tol) 
	         hopping(ih) = type20(ir,jr(ih),layer);
	      elseif(abs(d19-1)<tol)
	         hopping(ih) = type19(ir,jr(ih),layer);
	      elseif(abs(d21-1)<tol)
	         hopping(ih) = -type19(ir,jr(ih),layer);
	      elseif(abs(d19+1)<tol) 
	         hopping(ih) = -typem21(ir,jr(ih),layer);
	      elseif(abs(d21+1)<tol)
	         hopping(ih) = typem21(ir,jr(ih),layer);
	      else
	         hopping(ih) = 0.0;
	      end
           end
       end
    end
end     
end

