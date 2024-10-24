%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Routine to set TB parameters
% Adapted from matlab Shiang Fang matlab code:
% Ab initio tight-binding Hamiltonian for transition metal dichalcogenides
% by Shiang Fang, Rodrick Kuate Defo, Sharmila N. Shirodkar, Simon Lieu, Georgios A. Tritsaris, and Efthimios Kaxiras
% code version: July 2017
% Reference citation: Phys. Rev. B 92, 205108 (2015).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [onsite,type1,type2,type3,type4,type5,type7,type8,type10,type12,type14,...
     type16,typem16,type18,typem18,type19,type20,typem21,lsocM,lsocX] = set_parameters(tmdc,gw_par,sixth_nn)

% Read in parameters and classify them in terms of their type
TMDC_monolayer;
tbparms = mono_tbh_parm;
%vparms = mono_vint_parm;
%inter_par;

sq3=sqrt(3);
sq2=sqrt(2);

onsite = zeros(1,11);
type1  = zeros(11);
type2  = zeros(11);
type3  = zeros(11);
type4  = zeros(11);
type5  = zeros(11);
type7  = zeros(11);
type8  = zeros(11);
type10 = zeros(11);
type12 = zeros(11);
type14 = zeros(11);
type16 = zeros(11);
typem16= zeros(11);
type18 = zeros(11);
typem18= zeros(11);
type19 = zeros(11);
type20 = zeros(11);
typem21= zeros(11);
% 108 Independent parameters

gw1 = 0.0;
gw2 = gw1;
if(gw_par)
   gw1 = 0.3624;
   gw2 = -0.2512;
end

% ON-SITE
ep1=gw1 + tbparms(1);
ep3=gw2 + tbparms(2);
ep4=gw2 + tbparms(3);
ep6=gw1 + tbparms(4);
ep7=gw1 + tbparms(5);
ep9=gw2 + tbparms(6);
ep10=gw2 + tbparms(7);
ep2=ep1;
ep5=ep4;
ep8=ep7;
ep11=ep10;
onsite = [ep1,ep2,ep3,ep4,ep5,ep6,ep7,ep8,ep9,ep10,ep11];

gw1 = 1.0;
gw2 = gw1;
if(gw_par)
   gw1 = 1.4209;
   gw2 = 1.1738;
end

% TYPE 1 OFF-DIAGONAL
type1(3,5)  = gw2*tbparms(19);
type1(6,8)  = gw1*tbparms(20);
type1(9,11) = gw2*tbparms(21);
type1(1,2)  = gw1*tbparms(22);
type1(3,4)  = gw2*tbparms(23);
type1(4,5)  = gw2*tbparms(24);
type1(6,7)  = gw1*tbparms(25);
type1(7,8)  = gw1*tbparms(26);
type1(9,10) = gw2*tbparms(27);
type1(10,11)= gw2*tbparms(28);

type1 = type1 + type1';

% TYPE 1 DIAGONAL
type1(1,1)  = gw1*tbparms(8);
type1(2,2)  = gw1*tbparms(9);
type1(3,3)  = gw2*tbparms(10);
type1(4,4)  = gw2*tbparms(11);
type1(5,5)  = gw2*tbparms(12);
type1(6,6)  = gw1*tbparms(13);
type1(7,7)  = gw1*tbparms(14);
type1(8,8)  = gw1*tbparms(15);
type1(9,9)  = gw2*tbparms(16);
type1(10,10)= gw2*tbparms(17);
type1(11,11)= gw2*tbparms(18);

gw1 = 1.0;
if(gw_par)
   gw1 = 1.0773;
end
% TYPE 5 OFF-DIAGONAL ONLY
type5(4,1) = gw1*tbparms(29);
type5(3,2) = gw1*tbparms(30);
type5(5,2) = gw1*tbparms(31);
type5(9,6) = gw1*tbparms(32);
type5(11,6)= gw1*tbparms(33);
type5(10,7)= gw1*tbparms(34);
type5(9,8)=  gw1*tbparms(35);
type5(11,8)= gw1*tbparms(36);

gw1 = 1.0;
if(gw_par)
   gw1 = 1.1871;
end
if(sixth_nn)
  % TYPE 7 OFF-DIAGONAL ONLY
  type7(11,8)= gw1*tbparms(40); 
  type7(9,8)= gw1*tbparms(42);  
  type7(9,6) = gw1*tbparms(43); 
  type7(11,6)= gw1*tbparms(44); 
  type7(4,1) = gw1*tbparms(37); 
  type7(10,7)= gw1*tbparms(38); 
  type7(5,2) = gw1*tbparms(39); 
  type7(3,2) = gw1*tbparms(41); 
  
  % TYPE 10 OFF-DIAGONAL ONLY
  type10(4,1) = tbparms(45);
  type10(10,7)= tbparms(46);
  type10(5,2) = tbparms(47);
  type10(11,8)= tbparms(48);
  type10(3,2) = tbparms(49);
  type10(9,8)=  tbparms(50);
  type10(9,6) = tbparms(51);
  type10(11,6)= tbparms(52);
  type10(4,2) = tbparms(53);
  type10(10,8)= tbparms(54);
  type10(5,1) = tbparms(55);
  type10(11,7)= tbparms(56);
  type10(3,1) = tbparms(57);
  type10(9,7)=  tbparms(58);
  type10(10,6)= tbparms(59);
  
  type10 = type10+type10';
  % TYPE 18 OFF DIAGONAL
  type18(3,5)  = tbparms(71);
  type18(6,8)  = tbparms(72);
  type18(9,11)  = tbparms(73);
  
  type18 = type18 + type18';
  
  % TYPE 18 DIAGONAL
  type18(1,1)  = tbparms(60);
  type18(2,2)  = tbparms(61);
  type18(3,3)  = tbparms(62);
  type18(4,4)  = tbparms(63);
  type18(5,5)  = tbparms(64);
  type18(6,6)  = tbparms(65);
  type18(7,7)  = tbparms(66);
  type18(8,8)  = tbparms(67);
  type18(9,9)  = tbparms(68);
  type18(10,10)= tbparms(69);
  type18(11,11)= tbparms(70);
  
  % TYPE -18 OFF DIAGONAL
  typem18(3,5)  = tbparms(85);
  typem18(6,8)  = tbparms(86);
  typem18(9,11)  = tbparms(87);
  
  typem18 = typem18 + typem18';
  
  % TYPE -18 DIAGONAL
  typem18(1,1)  = tbparms(74);
  typem18(2,2)  = tbparms(75);
  typem18(3,3)  = tbparms(76);
  typem18(4,4)  = tbparms(77);
  typem18(5,5)  = tbparms(78);
  typem18(6,6)  = tbparms(79);
  typem18(7,7)  = tbparms(80);
  typem18(8,8)  = tbparms(81);
  typem18(9,9)  = tbparms(82);
  typem18(10,10)= tbparms(83);
  typem18(11,11)= tbparms(84);
  
  % TYPE 20 OFF DIAGONAL
  type20(3,5)  = tbparms(99);
  type20(6,8)  = tbparms(100);
  type20(9,11)  = tbparms(101);
  type20(1,2)  = tbparms(102);
  type20(3,4)  = tbparms(103);
  type20(4,5)  = tbparms(104);
  type20(6,7)  = tbparms(105);
  type20(7,8)  = tbparms(106);
  type20(9,10)  = tbparms(107);
  type20(10,11)= tbparms(108);
  
  type20 = type20 + type20';
  
  % TYPE 20 DIAGONAL
  type20(1,1)  = tbparms(88);
  type20(2,2)  = tbparms(89);
  type20(3,3)  = tbparms(90);
  type20(4,4)  = tbparms(91);
  type20(5,5)  = tbparms(92);
  type20(6,6)  = tbparms(93);
  type20(7,7)  = tbparms(94);
  type20(8,8)  = tbparms(95);
  type20(9,9)  = tbparms(96);
  type20(10,10)= tbparms(97);
  type20(11,11)= tbparms(98);
end  
% generate all remaining parameters

% TYPE 2 OFF-DIAGONAL

type2(3,5)=(sq3/2)*type1(3,4)-0.5*type1(3,5);
type2(6,8)=(sq3/2)*type1(6,7)-0.5*type1(6,8);
type2(9,11)=(sq3/2)*type1(9,10)-0.5*type1(9,11);

type2(1,2)=(sq3/4)*(type1(1,1)-type1(2,2))-type1(1,2);
type2(4,5)=(sq3/4)*(type1(4,4)-type1(5,5))-type1(4,5);
type2(7,8)=(sq3/4)*(type1(7,7)-type1(8,8))-type1(7,8);
type2(10,11)=(sq3/4)*(type1(10,10)-type1(11,11))-type1(10,11);

type2(3,4)=0.5*type1(3,4)+(sq3/2)*type1(3,5);
type2(6,7)=0.5*type1(6,7)+(sq3/2)*type1(6,8);
type2(9,10)=0.5*type1(9,10)+(sq3/2)*type1(9,11);

type2 = type2+type2';

% TYPE 2 DIAGONAL
type2(1,1)=0.25*type1(1,1)+0.75*type1(2,2);
type2(4,4)=0.25*type1(4,4)+0.75*type1(5,5);
type2(7,7)=0.25*type1(7,7)+0.75*type1(8,8);
type2(10,10)=0.25*type1(10,10)+0.75*type1(11,11);

type2(2,2)=0.75*type1(1,1)+0.25*type1(2,2);
type2(5,5)=0.75*type1(4,4)+0.25*type1(5,5);
type2(8,8)=0.75*type1(7,7)+0.25*type1(8,8);
type2(11,11)=0.75*type1(10,10)+0.25*type1(11,11);

type2(3,3)=type1(3,3);
type2(6,6)=type1(6,6);
type2(9,9)=type1(9,9);

% TYPE 3 OFF-DIAGONAL ONLY

type3(3,5)=-(sq3/2)*type1(3,4)-0.5*type1(3,5);
type3(6,8)=-(sq3/2)*type1(6,7)-0.5*type1(6,8);
type3(9,11)=-(sq3/2)*type1(9,10)-0.5*type1(9,11);

type3(1,2)=-(sq3/4)*(type1(1,1)-type1(2,2))-type1(1,2);
type3(4,5)=-(sq3/4)*(type1(4,4)-type1(5,5))-type1(4,5);
type3(7,8)=-(sq3/4)*(type1(7,7)-type1(8,8))-type1(7,8);
type3(10,11)=-(sq3/4)*(type1(10,10)-type1(11,11))-type1(10,11);

type3(3,4)=0.5*type1(3,4)-(sq3/2)*type1(3,5);
type3(6,7)=0.5*type1(6,7)-(sq3/2)*type1(6,8);
type3(9,10)=0.5*type1(9,10)-(sq3/2)*type1(9,11);

% TYPE 4 OFF-DIAGONAL ONLY
type4(4,1)=0.25*type5(4,1)+0.75*type5(5,2);
type4(10,7)=0.25*type5(10,7)+0.75*type5(11,8);

type4(5,2)=0.75*type5(4,1)+0.25*type5(5,2);
type4(11,8)=0.75*type5(10,7)+0.25*type5(11,8);

type4(5,1)=-(sq3/4)*type5(4,1)+(sq3/4)*type5(5,2);
type4(4,2)=type4(5,1);
type4(11,7)=-(sq3/4)*type5(10,7)+(sq3/4)*type5(11,8);
type4(10,8)=type4(11,7);

type4(3,1)=-(sq3/2)*type5(3,2);
type4(9,7)=-(sq3/2)*type5(9,8);

type4(3,2)=-0.5*type5(3,2);
type4(9,8)=-0.5*type5(9,8);

type4(9,6)=type5(9,6);
type4(10,6)=-sq3*type5(11,6)/2;
type4(11,6)=-type5(11,6)/2;


type3 = type3+type3';
type4 = type4+type4';
type5 = type5+type5';
type7 = type7+type7';

if(sixth_nn)
   alpha = [4, 7, 10];
   beta = [5, 8, 11];
   gamma = [3, 6, 9];
   
   for i = 1 : 3
       a = alpha(i); b = beta(i); c = gamma(i);
   
       typem16(a, a) = type18(a, a) / 4.0 + 3.0 * type18(b, b) / 4.0;
       typem16(b, b) = 3.0 * type18(a, a) / 4.0 + type18(b, b) / 4.0;
       typem16(c, c) = type18(c, c);
       typem16(c, a) = -sqrt(3) / 2.0 * type18(c, b);
       typem16(c, b) = -type18(c, b) / 2.0;
       typem16(a, b) = -sqrt(3) / 4.0 * (type18(a, a) - type18(b, b));
       
       typem16(b, a) = typem16(a, b); typem16(a, c) = typem16(c, a); typem16(b, c) = typem16(c, b);
       
       type16(a, a) = typem18(a, a) / 4.0 + 3.0 * typem18(b, b) / 4.0;
       type16(b, b) = 3.0 * typem18(a, a) / 4.0 + typem18(b, b) / 4.0;
       type16(c, c) = typem18(c, c);
       type16(c, a) = -sqrt(3) / 2.0 * typem18(c, b);
       type16(c, b) = -typem18(c, b) / 2.0;
       type16(a, b) = -sqrt(3) / 4.0 * (typem18(a, a) - typem18(b, b));
       
       type16(b, a) = type16(a, b); type16(a, c) = type16(c, a); type16(b, c) = type16(c, b);
       
       type19(a, a) = type20(a, a) / 4 + 3 * type20(b, b) / 4;
       type19(b, b) = 3 * type20(a, a) / 4 + type20(b, b) / 4;
       type19(c, c) = type20(c, c);
       type19(c, a) = -type20(c, a) / 2 - sqrt(3) / 2 * type20(c, b);
       type19(c, b) = sqrt(3) / 2 * type20(c, a) - type20(c, b) / 2;
       type19(a, b) = -sqrt(3) / 4 * (type20(a, a) - type20(b, b)) + type20(a, b);
   
       type19(b, a) = type19(a, b); type19(a, c) = type19(c, a); type19(b, c) = type19(c, b);
   
       typem19(c, a) = type20(c, a) / 2 - sqrt(3) / 2 * type20(c, b);
       typem19(c, b) = -sqrt(3) / 2 * type20(c, a) - type20(c, b) / 2;
       typem19(a, b) = -sqrt(3) / 4 * (type20(a, a) - type20(b, b)) - type20(a, b);
       typem19(b, a) = typem19(a, b); typem19(a, c) = typem19(c, a); typem19(b, c) = typem19(c, b);
   end
   
   a = 1; b = 2;
   typem16(a, a) = type18(a, a) / 4.0 + 3.0 * type18(b, b) / 4.0;
   typem16(b, b) = 3.0 * type18(a, a) / 4.0 + type18(b, b) / 4.0;
   typem16(a, b) = -sqrt(3) / 4.0 * (type18(a, a) - type18(b, b));
   typem16(b, a) = typem16(a, b);
   
   type16(a, a) = typem18(a, a) / 4.0 + 3.0 * typem18(b, b) / 4.0;
   type16(b, b) = 3.0 * typem18(a, a) / 4.0 + typem18(b, b) / 4.0;
   type16(a, b) = -sqrt(3) / 4.0 * (typem18(a, a) - typem18(b, b));
   type16(b, a) = type16(a, b);
   
   type19(a, a) = type20(a, a) / 4.0 + 3.0 * type20(b, b) / 4.0;
   type19(b, b) = 3 * type20(a, a) / 4 + type20(b, b) / 4;
   type19(a, b) = -sqrt(3) / 4.0 * (type20(a, a) - type20(b, b)) + type20(a, b);
   type19(b, a) = type19(a, b);
   
   typem19(a, b) = -sqrt(3) / 4.0 * (type20(a, a) - type20(b, b)) - type20(a, b);
   typem19(b, a) = typem19(a, b);
   
   alpha1 = [1, 7];
   beta1 = [2, 8];
   alpha2 = [4, 10];
   beta2 = [5, 11];
   gamma = [3, 9];
   
   for i = 1:2
       a1 = alpha1(i); a2 = alpha2(i); b1 = beta1(i); b2 = beta2(i); c2 = gamma(i);
   
       type8(a2, a1) = type7(a2, a1) / 4.0 + 3.0 * type7(b2, b1) / 4.0;
       type8(b2, b1) = 3.0 * type7(a2, a1) / 4.0 + type7(b2, b1) / 4.0;
       type8(b2, a1) = -sqrt(3) / 4.0 * (type7(a2, a1)-type7(b2, b1));
       type8(a2, b1) = type8(b2, a1);
       type8(c2, a1) = -sqrt(3) / 2.0 * type7(c2, b1);
       type8(c2, b1) = -type7(c2, b1) / 2.0;
       
       type12(a2, a1) = type10(a2, a1) / 4.0 + 3.0 * type10(b2, b1) / 4.0 - sqrt(3) / 4.0 * (type10(a2,b1) + type10(b2, a1));
       type12(b2, b1) = 3.0 * type10(a2, a1) / 4.0 + type10(b2, b1) / 4.0 + sqrt(3) / 4.0 * (type10(a2,b1) + type10(b2, a1));
       type12(b2, a1) = sqrt(3) / 4 * (type10(a2, a1) - type10(b2, b1)) + type10(b2, a1) / 4 - 3 * type10(a2, b1) / 4;
       type12(a2, b1) = sqrt(3) / 4 * (type10(a2, a1) - type10(b2, b1)) - 3 * type10(b2, a1) / 4 + type10(a2, b1) / 4;
       type12(c2, a1) = -type10(c2, a1) / 2 + sqrt(3) / 2 * type10(c2, b1);
       type12(c2, b1) = -sqrt(3) / 2 * type10(c2, a1) - type10(c2, b1) / 2;
   
       type14(a2, a1) = type10(a2, a1) / 4.0 + 3.0 * type10(b2, b1) / 4.0 + sqrt(3) / 4.0 * (type10(a2,b1) + type10(b2, a1));
       type14(b2, b1) = 3.0 * type10(a2, a1) / 4.0 + type10(b2, b1) / 4.0 - sqrt(3) / 4.0 * (type10(a2,b1) + type10(b2, a1));
       type14(b2, a1) = -sqrt(3) / 4 * (type10(a2, a1) - type10(b2, b1)) + type10(b2, a1) / 4 - 3 * type10(a2, b1) / 4;
       type14(a2, b1) = -sqrt(3) / 4 * (type10(a2, a1) - type10(b2, b1)) - 3 * type10(b2, a1) / 4 + type10(a2, b1) / 4;
       type14(c2, a1) = -type10(c2, a1) / 2 - sqrt(3) / 2 * type10(c2, b1);
       type14(c2, b1) = sqrt(3) / 2 * type10(c2, a1) - type10(c2, b1) / 2;
   end
   
   type8(9, 6) = type7(9, 6);
   type8(10, 6) = -sqrt(3) / 2.0 * type7(11, 6);
   type8(11, 6) = -type7(11, 6) / 2.0;
   
   type12(9, 6) = type10(9, 6);
   type12(10, 6) = -type10(10, 6) / 2 + sqrt(3) * type10(11, 6) / 2;
   type12(11, 6) = -sqrt(3) / 2 * type10(10, 6) - type10(11, 6) / 2;
   
   type14(9, 6) = type10(9, 6);
   type14(10, 6) = -type10(10, 6) / 2 - sqrt(3) * type10(11, 6) / 2;
   type14(11, 6) = sqrt(3) / 2 * type10(10, 6) - type10(11, 6) / 2;
   
   type8 = type8 + type8';
   type12 = type12 + type12';
   type14 = type14 + type14';
   
   % forgotten non zero component I am such an idiot
   even1 = [3, 6, 9]; even2 = [5, 8, 11];
   odd1 = [1, 3, 4, 6, 7, 9, 10]; odd2 = [2, 4, 5, 7, 8, 10, 11];
   for i = 1:3
       typem21(even1(i), even2(i)) = typem19(even1(i), even2(i));
   end
   for i = 1:7
       typem21(odd1(i), odd2(i)) = -typem19(odd1(i), odd2(i));%-typem19(odd1(i), odd2(i));
   end
   
   typem21 = typem21 + typem21';
else
  % TYPE 7 OFF-DIAGONAL ONLY
  type7(9,6) = gw1*tbparms(37); 
  type7(11,6)= gw1*tbparms(38);
  type7(9,8) =  gw1*tbparms(39);  
  type7(11,8)= gw1*tbparms(40); 
end 
%Vsigma=vparms(1);
%Vpi=vparms(2);
%R=vparms(3:4);
%eta=vparms(5:6);

lsocM = soc_m;
lsocX = soc_x;

end

