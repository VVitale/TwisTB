%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slater-Koster parameters for p-p and p-d
% hoppings calculated at different interlayer
% seperations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% multilayer  - Boolean that gives whether calc is 
%            a multilayer calculation
% tmdc     - 2x1 array that has indices of top 
%            and bottom TMD, i.e. 1-MoS2,4-WSe2
%
% Outputs:
% pp_vint_z     - z-distance between p-p orbitals
%           that the parameters are calculated for              
% pp_vint_parm  - p-p SK parameters at different 
%               interlayer seperations
% pd_vint_lay1_z    - z-distance between bottom d 
%                   and top p orbital
% pd_vint_lay1_parm - p-d SK for bottom d and top p
% pd_vint_lay2_z    - z-dist from top d to bottom p
% pd_vint_lay2_parm - p-d SK for top d and bottom p
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Kemal Atalar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pp_vint_z,pp_vint_parm,pd_vint_lay1_z,pd_vint_lay1_parm,pd_vint_lay2_z,pd_vint_lay2_parm] = inter_par(multilayer,nlayers,tmdc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Homobilayer p-p and p-d parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MoS2-MoS2
% MoS2 p-p interlayer parameters, at distances zpp
% p-p
MoS2_MoS2_zpp={2.6170,2.7170,2.8170,2.9170,3.0170,3.1170,3.2170,3.3170,3.4170};
MoS2_MoS2_pp_inter={
      [4.07660,-1.27092,2.82371,2.61189,3.54713,5.41905];
      [3.72336,-1.11170,2.88775,2.68782,3.59004,5.49683];
      [3.39650,-0.98342,2.95225,2.75873,3.62986,5.52781];
      [3.09684,-0.88484,3.01667,2.82178,3.66561,5.49774];
      [2.82420,-0.81786,3.08049,2.87201,3.69634,5.38604];
      [2.57738,-0.79236,3.14336,2.89917,3.72152,5.15810];
      [2.35482,-0.84288,3.20489,2.87939,3.74068,4.75358];
      [2.15466,-1.13777,3.26482,2.73898,3.75364,4.04510];
      [1.97483,-6.14163,3.32293,1.95866,3.76049,2.50090];
};
% p-d
MoS2_MoS2_lay1_zpd={4.2085,4.3085,4.4085,4.5085,4.6085,4.7085,4.8085,4.9085,5.0085};
MoS2_MoS2_pd_lay1_inter={
      [-603.92946,-1779650.51956,-5.84309,-10.33203,2.38593,1.71046,-4.78384,-2.38800];
      [-343.84624,-1570234.52372,-5.45520,-10.18617,2.24456,1.68921,-4.24656,-2.41207];
      [-364.96329,-1395721.13766,-5.46461,-10.05024,2.17058,1.66622,-4.04457,-2.42325];
      [-464.61193,-1252974.52153,-5.58757,-9.92533,2.11978,1.64161,-3.95931,-2.42197];
      [-641.97215,-1139297.26973,-5.76010,-9.81268,2.08069,1.61569,-3.93348,-2.40953];
      [-925.68538,-1052163.41409,-5.95716,-9.71349,2.04866,1.58878,-3.94345,-2.38747];
      [-1367.89540,-986090.36204,-6.16730,-9.62684,2.02163,1.56115,-3.97876,-2.35689];
      [-2043.22680,-933160.09895,-6.38203,-9.54941,1.99763,1.53287,-4.02876,-2.31791];
      [-3068.96839,-893373.37926,-6.59821,-9.48202,1.97636,1.50429,-4.09214,-2.27231];
};
% p-d from top metal to bottom chalcogen
MoS2_MoS2_lay2_zpd={4.2085,4.3085,4.4085,4.5085,4.6085,4.7085,4.8085,4.9085,5.0085};
%MoS2_MoS2_pd_lay2_inter={
%      [12047.22193,289760.18349,-7.76203,-9.22269,2.49255,1.73230,-5.44840,-2.34302];
%      [5005.45789,259321.06262,-7.15810,-9.09696,2.36120,1.69373,-4.97630,-2.28485];
%      [3541.24618,239809.18166,-6.90062,-8.99556,2.25864,1.65738,-4.63211,-2.23112];
%      [3424.29031,231435.31570,-6.84214,-8.92367,2.17857,1.62342,-4.39311,-2.18428];
%      [3873.84573,234498.87118,-6.88169,-8.88395,2.11262,1.59217,-4.21992,-2.14751];
%      [4776.77765,250187.01687,-6.97384,-8.87724,2.05641,1.56392,-4.09182,-2.12363];
%      [6180.75552,281129.48074,-7.09456,-8.90281,2.00728,1.53886,-3.99588,-2.11487];
%      [8220.24204,331595.04889,-7.23099,-8.95792,1.96359,1.51699,-3.92416,-2.12214];
%      [11106.77876,406128.90585,-7.37599,-9.03562,1.92429,1.49780,-3.87166,-2.14297];
%};
% Symmetrize so same as layer1 just different sign
MoS2_MoS2_pd_lay2_inter={
      [603.92946,1779650.51956,-5.84309,-10.33203,2.38593,1.71046,-4.78384,-2.38800];
      [343.84624,1570234.52372,-5.45520,-10.18617,2.24456,1.68921,-4.24656,-2.41207];
      [364.96329,1395721.13766,-5.46461,-10.05024,2.17058,1.66622,-4.04457,-2.42325];
      [464.61193,1252974.52153,-5.58757,-9.92533,2.11978,1.64161,-3.95931,-2.42197];
      [641.97215,1139297.26973,-5.76010,-9.81268,2.08069,1.61569,-3.93348,-2.40953];
      [925.68538,1052163.41409,-5.95716,-9.71349,2.04866,1.58878,-3.94345,-2.38747];
      [1367.89540,986090.36204,-6.16730,-9.62684,2.02163,1.56115,-3.97876,-2.35689];
      [2043.22680,933160.09895,-6.38203,-9.54941,1.99763,1.53287,-4.02876,-2.31791];
      [3068.96839,893373.37926,-6.59821,-9.48202,1.97636,1.50429,-4.09214,-2.27231];
};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MoSe2-MoSe2
% MoSe2 p-p interlayer parameters, at distance zpp
MoSe2_MoSe2_zpp={2.4850,2.5850,2.6850,2.7850,2.8850,2.9850,3.0850,3.1850,3.2850};
MoSe2_MoSe2_pp_inter={
      [4.64287,-1.97036,2.88111,2.51311,3.64215,5.04578];
      [4.30810,-1.74472,2.93987,2.58816,3.68439,5.15746];
      [3.98402,-1.55447,3.00067,2.66023,3.72644,5.23449];
      [3.67379,-1.39773,3.06319,2.72775,3.76824,5.26873];
      [3.38136,-1.27678,3.12674,2.78749,3.80830,5.24518];
      [3.10875,-1.19704,3.19075,2.83428,3.84547,5.14403];
      [2.85655,-1.17622,3.25482,2.85730,3.87906,4.93086];
      [2.62487,-1.27587,3.31852,2.82962,3.90822,4.53871];
      [2.41308,-1.84610,3.38148,2.65430,3.93244,3.79729];
};
% p-d
MoSe2_MoSe2_lay1_zpd={4.1425,4.2425,4.3425,4.4425,4.5425,4.6425,4.7425,4.8425,4.9425};
MoSe2_MoSe2_pd_lay1_inter={
      [-256511.97369,-3029864.72589,-9.65586,-10.55821,2.77041,1.60534,-6.69552,-1.92380];
      [-963.79409,-2558694.29419,-6.01469,-10.37853,2.34040,1.60144,-4.76841,-2.01520];
      [-261.08635,-2124434.50291,-5.15944,-10.19463,2.14096,1.59434,-3.89911,-2.08968];
      [-258.77673,-1758627.66917,-5.12426,-10.01479,2.06770,1.58340,-3.67542,-2.14376];
      [-321.15559,-1461530.02272,-5.22879,-9.84279,2.02166,1.56888,-3.59444,-2.17784];
      [-435.47841,-1227167.34302,-5.38691,-9.68203,1.98747,1.55123,-3.57682,-2.19337];
      [-617.97738,-1046383.26612,-5.57142,-9.53509,1.96047,1.53094,-3.59807,-2.19235];
      [-898.07787,-906992.40603,-5.76877,-9.40211,1.93785,1.50849,-3.64323,-2.17653];
      [-1323.34179,-799324.49076,-5.97279,-9.28276,1.91843,1.48430,-3.70583,-2.14770];
};
MoSe2_MoSe2_lay2_zpd={4.1425,4.2425,4.3425,4.4425,4.5425,4.6425,4.7425,4.8425,4.9425};
%MoSe2_MoSe2_pd_lay2_inter={
%      [54741.17396,473505.00675,-8.64776,-9.42383,2.55592,1.67580,-5.76544,-2.10494];
%      [18015.87686,393579.72937,-7.87730,-9.24521,2.43864,1.64553,-5.37631,-2.07693];
%      [5930.54257,328513.57203,-7.12697,-9.07502,2.30607,1.61651,-4.87872,-2.04830];
%      [3452.13429,279399.73736,-6.74722,-8.92138,2.19787,1.58794,-4.48527,-2.01673];
%      [3040.89769,246368.16531,-6.62958,-8.79392,2.11631,1.55999,-4.22017,-1.98449];
%      [3301.69464,228384.80034,-6.64179,-8.69998,2.05167,1.53329,-4.03798,-1.95614];
%      [3999.33108,224080.44434,-6.72054,-8.64257,1.99793,1.50850,-3.90911,-1.93633];
%      [5127.27344,232876.92742,-6.83289,-8.62116,1.95133,1.48608,-3.81449,-1.92851];
%      [6816.39297,256011.47478,-6.96610,-8.63408,1.91061,1.46638,-3.74821,-1.93543];
%};
MoSe2_MoSe2_pd_lay2_inter={
      [256511.97369,3029864.72589,-9.65586,-10.55821,2.77041,1.60534,-6.69552,-1.92380];
      [963.79409,2558694.29419,-6.01469,-10.37853,2.34040,1.60144,-4.76841,-2.01520];
      [261.08635,2124434.50291,-5.15944,-10.19463,2.14096,1.59434,-3.89911,-2.08968];
      [258.77673,1758627.66917,-5.12426,-10.01479,2.06770,1.58340,-3.67542,-2.14376];
      [321.15559,1461530.02272,-5.22879,-9.84279,2.02166,1.56888,-3.59444,-2.17784];
      [435.47841,1227167.34302,-5.38691,-9.68203,1.98747,1.55123,-3.57682,-2.19337];
      [617.97738,1046383.26612,-5.57142,-9.53509,1.96047,1.53094,-3.59807,-2.19235];
      [898.07787,906992.40603,-5.76877,-9.40211,1.93785,1.50849,-3.64323,-2.17653];
      [1323.34179,799324.49076,-5.97279,-9.28276,1.91843,1.48430,-3.70583,-2.14770];
};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WS2-WS2
% p-p
WS2_WS2_zpp={0.0, 2.6180,2.7180,2.8180,2.9180,3.0180,3.1180,3.2180,3.3180};
WS2_WS2_pp_inter={
      [0.00000, 0.00000,0.00000,0.00000,0.00000,0.00000];
      [3.90795,-1.02097,2.84655,2.70446,3.54628,6.13100];
      [3.55262,-0.85858,2.91511,2.79938,3.59832,6.34796];
      [3.25312,-0.83053,2.97958,2.83687,3.63551,6.09974];
      [2.98662,-0.71419,3.03971,2.92167,3.65795,6.18028];
      [2.71910,-0.70953,3.10544,2.94457,3.69106,5.81570];
      [2.46881,-0.70544,3.17253,2.96179,3.72388,5.45537];
      [2.24472,-0.80992,3.23794,2.90569,3.74936,4.81419];
      [2.04267,-1.53765,3.30210,2.59232,3.76926,3.61969];
};
% p-d
WS2_WS2_lay1_zpd={0.0,4.2090,4.3090,4.4090,4.5090,4.6090,4.7090,4.8090,4.9090};
WS2_WS2_pd_lay1_inter={
      [  0.00000, 0.0000000000, 0.00000, 0.00000,0.00000,0.00000, 0.00000, 0.00000];
      [-23.10206,-776085.48403,-3.74376,-9.76056,2.15287,1.64544,-3.58399,-2.12767];
      [-32.19062,-458451.42386,-3.92755,-9.36903,2.11528,1.63843,-3.55945,-2.19406];
      [-33.65918,-470319.04386,-3.92263,-9.32423,2.05210,1.62964,-3.39895,-2.26762];
      [-54.15999,-451722.17442,-4.19778,-9.24098,2.04484,1.61984,-3.53173,-2.32859];
      [-82.52882,-323330.60834,-4.43524,-8.98965,2.02867,1.59860,-3.62570,-2.32108];
      [-130.27112,-270118.32503,-4.69509,-8.83595,2.01870,1.57545,-3.75169,-2.30628];
      [-206.17122,-230800.18120,-4.95463,-8.70053,2.01033,1.54952,-3.88818,-2.27359];
      [-325.19473,-206099.04666,-5.21010,-8.59527,2.00271,1.52291,-4.02980,-2.23436];
};
WS2_WS2_lay2_zpd={0.0,4.2090,4.3090,4.4090,4.5090,4.6090,4.7090,4.8090,4.9090};
%WS2_WS2_pd_lay2_inter={
%      [209.88183,167644.23550,-5.15807,-8.84587,2.29525,1.72570,-4.44468,-2.33938];
%      [226.46919,103125.61263,-5.17695,-8.48554,2.23549,1.67662,-4.30026,-2.21020];
%      [256.53420,123326.76811,-5.21567,-8.53826,2.15738,1.66345,-4.07366,-2.27341];
%      [312.79674,90606.41496,-5.30613,-8.30205,2.09824,1.62092,-3.93162,-2.16719];
%      [482.50878,96103.63753,-5.54405,-8.28819,2.06820,1.60094,-3.94340,-2.18278];
%      [696.91897,90852.47848,-5.73932,-8.21119,2.03188,1.57212,-3.91951,-2.14602];
%      [1016.98786,95791.24350,-5.93930,-8.20157,1.99980,1.54928,-3.91671,-2.14287];
%      [1485.62512,107279.99100,-6.13878,-8.22967,1.97074,1.52912,-3.92745,-2.15477];
%};
WS2_WS2_pd_lay2_inter={
      [0.00000, 0.0000000000, 0.00000, 0.00000,0.00000,0.00000, 0.00000, 0.00000];
      [23.10206,776085.48403,-3.74376,-9.76056,2.15287,1.64544,-3.58399,-2.12767];
      [32.19062,458451.42386,-3.92755,-9.36903,2.11528,1.63843,-3.55945,-2.19406];
      [33.65918,470319.04386,-3.92263,-9.32423,2.05210,1.62964,-3.39895,-2.26762];
      [54.15999,451722.17442,-4.19778,-9.24098,2.04484,1.61984,-3.53173,-2.32859];
      [82.52882,323330.60834,-4.43524,-8.98965,2.02867,1.59860,-3.62570,-2.32108];
      [130.27112,270118.32503,-4.69509,-8.83595,2.01870,1.57545,-3.75169,-2.30628];
      [206.17122,230800.18120,-4.95463,-8.70053,2.01033,1.54952,-3.88818,-2.27359];
      [325.19473,206099.04666,-5.21010,-8.59527,2.00271,1.52291,-4.02980,-2.23436];
};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WSe2-WSe2
% p-p
WSe2_WSe2_zpp={2.4855,2.5855,2.6855,2.7855,2.8855,2.9855,3.0855,3.1855,3.2855};
WSe2_WSe2_pp_inter={
      [3.77532,-1.04153,3.02031,2.75682,3.94867,7.53891];
      [4.26945,-1.01179,2.92522,2.79872,3.55666,7.16000];
      [3.86290,-0.83403,3.00973,2.91045,3.66829,7.48761];
      [3.80744,-0.77887,3.02406,2.96963,3.59039,7.27066];
      [3.47484,-1.01937,3.09585,2.89059,3.64621,5.94986];
      [2.94072,-0.83752,3.22501,2.99760,3.83990,6.14819];
      [2.69241,-0.73659,3.29391,3.07555,3.88085,6.10354];
      [2.44160,-0.70889,3.36929,3.10991,3.93751,5.76339];
      [2.24240,-0.71238,-3.43370,3.11388,3.96233,5.27715];
};
% p-d
WSe2_WSe2_lay1_zpd={4.1427,4.2428,4.3428,4.4428,4.5427,4.6427,4.7428,4.8428,4.9428};
WSe2_WSe2_pd_lay1_inter={
      [-1.43942,-32854.36330,-2.04875,-7.74319,1.85235,1.46603,-2.08664,-1.34319];
      [-10789.83475,-392290.36323,-7.56037,-9.16069,2.65419,1.58554,-6.40885,-1.98108];
      [-5.60114,-49047.07589,-2.78930,-7.86671,1.90634,1.46231,-2.67179,-1.47574];
      [-31.11607,-108067.41025,-3.82395,-8.27466,1.98398,1.52807,-3.24441,-1.88891];
      [-55.07976,-88047.81379,-4.15009,-8.10147,1.98283,1.52002,-3.39476,-1.93241];
      [-59.08045,-130047.22501,-4.15073,-8.28552,1.98120,1.52868,-3.55170,-2.08566];
      [-132.78847,-82028.55823,-4.62079,-7.97348,2.00286,1.50138,-3.84855,-2.00904];
      [-204.02836,-63134.21420,-4.85034,-7.78223,1.99551,1.47923,-3.98013,-1.98490];
      [-302.88276,-84612.37927,-5.06082,-7.91268,1.99141,1.47669,-4.13219,-2.08904];
};
WSe2_WSe2_lay2_zpd={4.1427,4.2428,4.3428,4.4428,4.5427,4.6427,4.7428,4.8428,4.9428};
%WSe2_WSe2_pd_lay2_inter={
%      [192565.67948,72008.03493,-9.45778,-8.22398,2.70831,1.58655,-6.56236,-1.71823];
%      [22.34529,2323.81797,-3.70886,-6.11380,1.96279,1.41167,-2.99557,-0.86481];
%      [19961.94716,18095.12917,-7.87704,-7.27158,2.47010,1.55980,-5.78525,-1.74376];
%      [417.85268,10055.63430,-5.44269,-6.88187,2.12659,1.49517,-4.14112,-1.49584];
%      [458.38573,8204.45678,-5.45630,-6.72209,2.07784,1.46237,-4.03679,-1.40530];
%      [943.68324,5008.61352,-5.86091,-6.40276,2.07823,1.40765,-4.20715,-1.17885];
%      [799.08999,5540.77673,-5.71816,-6.42832,2.00493,1.36750,-3.96812,-1.04320];
%      [1240.21518,6528.55821,-5.94685,-6.49268,1.98547,1.35058,-4.02690,-1.04420];
%      [1316.28407,8431.93548,-5.94805,-6.61164,1.93823,1.33810,-3.92063,-1.07909];
%};
WSe2_WSe2_pd_lay2_inter={
      [1.43942,32854.36330,-2.04875,-7.74319,1.85235,1.46603,-2.08664,-1.34319];
      [10789.83475,392290.36323,-7.56037,-9.16069,2.65419,1.58554,-6.40885,-1.98108];
      [5.60114,49047.07589,-2.78930,-7.86671,1.90634,1.46231,-2.67179,-1.47574];
      [31.11607,108067.41025,-3.82395,-8.27466,1.98398,1.52807,-3.24441,-1.88891];
      [55.07976,88047.81379,-4.15009,-8.10147,1.98283,1.52002,-3.39476,-1.93241];
      [59.08045,130047.22501,-4.15073,-8.28552,1.98120,1.52868,-3.55170,-2.08566];
      [132.78847,82028.55823,-4.62079,-7.97348,2.00286,1.50138,-3.84855,-2.00904];
      [204.02836,63134.21420,-4.85034,-7.78223,1.99551,1.47923,-3.98013,-1.98490];
      [302.88276,84612.37927,-5.06082,-7.91268,1.99141,1.47669,-4.13219,-2.08904];
};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Heterobilayer p-p and p-d parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MoS2-MoSe2
%MoS2-MoSe2 p-p interlayer parameters for S-Se, at distance zpp
MoS2_MoSe2_zpp={2.6170,2.7170,2.8170,2.9170,3.0170,3.1170,3.2170,3.3170,3.4170};
MoS2_MoSe2_pp_inter={
      [3.98374,-1.39320,2.93078,2.64108,3.83778,5.78948];
      [3.65139,-1.22388,2.99599,2.71836,3.88590,5.88571];
      [3.33959,-1.08019,3.06245,2.79356,3.93278,5.94871];
      [3.04770,-0.96824,3.13013,2.86171,3.97854,5.94317];
      [2.78078,-0.88325,3.19769,2.92202,4.01989,5.87291];
      [2.53809,-0.82947,3.26472,2.96877,4.05585,5.71308];
      [2.31797,-0.82565,3.33103,2.98767,4.08605,5.40888];
      [2.11763,-0.92931,3.39678,2.94595,4.11172,4.87775];
      [1.93624,-1.57638,3.46157,2.70227,4.13200,3.87194];
};
% p-d
MoS2_MoSe2_lay1_zpd={4.2085,4.3085,4.4085,4.5085,4.6085,4.7085,4.8085,4.9085,5.0085};
MoS2_MoSe2_pd_lay1_inter={
      [-481.46531,-1404835.68798,-5.64875,-10.11822,2.39656,1.72040,-4.86590,-2.45147];
      [-269.02491,-1129372.26459,-5.24983,-9.91557,2.25680,1.70745,-4.33294,-2.51033];
      [-291.65874,-917343.86368,-5.27178,-9.72520,2.19021,1.69494,-4.16877,-2.56811];
      [-374.00964,-791175.14822,-5.39867,-9.57762,2.14425,1.67530,-4.10870,-2.59099];
      [-515.54462,-708306.84342,-5.56932,-9.45760,2.10950,1.65383,-4.10633,-2.60212];
      [-740.27389,-656104.72801,-5.76359,-9.36222,2.08116,1.62986,-4.13701,-2.59783];
      [-1090.58909,-619804.33089,-5.97163,-9.28209,2.05760,1.60471,-4.19265,-2.58375];
      [-1622.76546,-599343.49049,-6.18379,-9.21881,2.03671,1.57911,-4.26186,-2.56353];
      [-2428.01849,-587000.64732,-6.39737,-9.16584,2.01815,1.55300,-4.34281,-2.53647];
};
MoS2_MoSe2_lay2_zpd={4.2085,4.3085,4.4085,4.5085,4.6085,4.7085,4.8085,4.9085,5.0085};
MoS2_MoSe2_pd_lay2_inter={
      [3087.22313,209660.62008,-6.76154,-8.92619,2.26907,1.66230,-4.23761,-1.90574];
      [2726.47013,181872.32958,-6.63639,-8.77858,2.18041,1.62549,-3.95979,-1.83942];
      [2923.51876,160786.48307,-6.63671,-8.64748,2.10928,1.58803,-3.76299,-1.76529];
      [3619.75490,157343.56713,-6.72701,-8.58145,2.05435,1.55776,-3.64510,-1.72546];
      [4667.29089,161288.29647,-6.84210,-8.54680,2.00507,1.52751,-3.55075,-1.68406];
      [6244.43030,176521.98458,-6.97897,-8.55341,1.96213,1.50108,-3.48550,-1.66197];
      [8557.99722,206025.42144,-7.12944,-8.59938,1.92438,1.47917,-3.44430,-1.66432];
      [11880.25728,252047.76247,-7.28667,-8.67380,1.89068,1.46099,-3.42176,-1.68736];
      [16633.50221,321960.83399,-7.44811,-8.77372,1.86044,1.44650,-3.41518,-1.73164];
};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MoS2-WS2
% p-p
MoS2_WS2_zpp={0.0,2.6170,2.7170,2.8170,2.9170,3.0170,3.1170,3.2170,3.3170};
MoS2_WS2_pp_inter={
      [0.00000, 0.00000,0.00000,0.00000,0.00000,0.00000];
      [4.06380,-1.18176,2.82369,2.64417,3.52230,5.67773];
      [3.70887,-1.03778,2.88870,2.71923,3.56647,5.74120];
      [3.38690,-0.92649,2.95272,2.78696,3.60350,5.73562];
      [3.08863,-0.84453,3.01715,2.84485,3.63759,5.65153];
      [2.81252,-0.79403,3.08213,2.88763,3.66961,5.47194];
      [2.56293,-0.79742,3.14598,2.89765,3.69566,5.13566];
      [2.33596,-0.91543,3.20914,2.83965,3.71731,4.57475];
      [2.13429,-1.58087,3.26988,2.57589,3.73137,3.60418];
};
% p-d
MoS2_WS2_lay1_zpd={0.0,4.2085,4.3085,4.4085,4.5085,4.6085,4.7085,4.8085,4.9085};
MoS2_WS2_pd_lay1_inter={
      [   0.00000, 0.0000000000, 0.00000, 0.00000,0.00000,0.00000, 0.00000, 0.00000];
      [-101.53749,-1249557.66515,-4.70353,-10.10930,2.20441,1.66213,-3.83375,-2.17052];
      [-125.41332,-1065320.06172,-4.81196,-9.94367,2.14758,1.65081,-3.71925,-2.23078];
      [-176.75798,-900530.27937,-5.00270,-9.77834,2.11106,1.63770,-3.71020,-2.28011];
      [-260.88135,-794977.68169,-5.22072,-9.64535,2.08217,1.62038,-3.74040,-2.30773];
      [-393.48140,-722044.32921,-5.45037,-9.53419,2.05839,1.60147,-3.79768,-2.32571];
      [-599.20020,-686491.75018,-5.68399,-9.45443,2.03901,1.58258,-3.87848,-2.34273];
      [-916.87016,-657386.77119,-5.91873,-9.38321,2.02061,1.55887,-3.96310,-2.33182];
      [-1377.08927,-604772.07253,-6.14066,-9.29144,2.00345,1.53201,-4.05461,-2.30085];
};
MoS2_WS2_lay2_zpd={0.0,4.2085,4.3085,4.4085,4.5085,4.6085,4.7085,4.8085,4.9085};
MoS2_WS2_pd_lay2_inter={
      [0.00000, 0.0000000000, 0.00000, 0.00000,0.00000,0.00000, 0.00000, 0.00000];
      [1341.56538,238309.48499,-6.32186,-9.05706,2.41605,1.75465,-5.07957,-2.47022];
      [728.80699,213731.94126,-5.89902,-8.92854,2.28997,1.71886,-4.61274,-2.42571];
      [711.54404,159172.67650,-5.84876,-8.69067,2.20868,1.68174,-4.36587,-2.35695];
      [862.47829,148943.10501,-5.93414,-8.59798,2.14709,1.65107,-4.21627,-2.32419];
      [1124.74577,135347.16420,-6.06441,-8.49247,2.09450,1.61610,-4.10858,-2.26222];
      [1542.73285,142332.37611,-6.22441,-8.47702,2.04870,1.59051,-4.03423,-2.25126];
      [2185.93680,158291.61947,-6.40278,-8.49689,2.00900,1.56725,-3.98846,-2.25125];
      [3056.84660,189686.14961,-6.57186,-8.56108,1.97777,1.55455,-3.98757,-2.30909];
};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MoS2-WSe2
% p-p
MoS2_WSe2_zpp={2.6180,2.7180,2.8180,2.9180,3.0180,3.1180,3.2180,3.3180};
MoS2_WSe2_pp_inter={
      [3.91409,-1.35682,2.93643,2.64582,3.79269,5.78831];
      [3.57585,-1.19066,3.00534,2.72403,3.84974,5.87604];
      [3.23939,-1.06245,3.08019,2.79461,3.91958,5.88730];
      [2.97585,-0.97945,3.14392,2.85207,3.95293,5.80485];
      [2.71780,-0.93666,3.21223,2.89085,4.00005,5.59080];
      [2.48808,-0.92244,3.27770,2.91545,4.03170,5.32233];
      [2.26266,-1.02033,3.34786,2.88116,4.07324,4.81768];
      [2.05294,-1.59811,3.41947,2.67207,4.11589,3.89001];
};
% p-d
MoS2_WSe2_lay1_zpd={4.2090,4.3090,4.4090,4.5090,4.6090,4.7090,4.8090,4.9090};
MoS2_WSe2_pd_lay1_inter={
      [-123.26003,-462735.34292,-4.79111,-9.41396,2.28193,1.69010,-4.28515,-2.31516];
      [-139.48539,-477529.62972,-4.84058,-9.37465,2.22204,1.68738,-4.15920,-2.41737];
      [-196.08786,-480521.79320,-5.02784,-9.32149,2.18389,1.68567,-4.14412,-2.52689];
      [-276.24420,-406269.30687,-5.21569,-9.16491,2.15228,1.67095,-4.16440,-2.57120];
      [-390.05640,-387257.35931,-5.40096,-9.08741,2.12314,1.65860,-4.19414,-2.62587];
      [-567.38646,-394815.91925,-5.60540,-9.05277,2.10137,1.63462,-4.26420,-2.62669];
      [-838.56430,-386937.56173,-5.81688,-8.99705,2.08150,1.61597,-4.34364,-2.64801];
      [-1237.57852,-402347.13683,-6.02577,-8.97836,2.06685,1.60028,-4.45227,-2.68534];
};
MoS2_WSe2_lay2_zpd={4.2090,4.3090,4.4090,4.5090,4.6090,4.7090,4.8090,4.9090};
MoS2_WSe2_pd_lay2_inter={
      [78.21324,75277.72675,-4.37382,-8.26552,2.04440,1.60989,-3.09675,-1.65551];
      [120.43159,22725.58446,-4.61523,-7.47206,2.01462,1.49640,-3.07965,-1.15616];
      [289.55946,25943.39302,-5.15163,-7.49687,1.99945,1.45602,-3.15810,-1.03911];
      [490.08384,67503.60243,-5.44378,-8.02239,1.97288,1.52379,-3.17727,-1.55594];
      [769.06467,64146.69064,-5.68538,-7.94610,1.93025,1.49415,-3.12333,-1.50253];
      [1040.90717,70251.06951,-5.83011,-7.95330,1.91021,1.47264,-3.16774,-1.50069];
      [1647.66756,97920.40577,-6.07426,-8.10364,1.89164,1.47501,-3.22716,-1.64168];
      [2596.39086,64394.36092,-6.31371,-7.81547,1.87996,1.39896,-3.32659,-1.28851];
};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WS2-MoSe2
% p-p
WS2_MoSe2_zpp={2.6180,2.7180,2.8180,2.9180,3.0180,3.1180,3.2180,3.3180,3.4180};
WS2_MoSe2_pp_inter={
      [3.91791,-1.03906,2.93822,2.75304,3.82989,6.80133];
      [3.50152,-0.96067,3.02366,2.81212,3.93860,6.68807];
      [3.20957,-0.82222,3.09033,2.90183,3.98619,6.85931];
      [2.92665,-0.72550,3.16029,2.97926,4.03748,6.89229];
      [2.66852,-0.65330,3.22964,3.04756,4.08174,6.80656];
      [2.43023,-0.61616,3.29922,3.09516,4.12223,6.53321];
      [2.20681,-0.60285,3.37025,3.12354,4.16368,6.14019];
      [2.00326,-0.63959,3.44085,3.10983,4.20112,5.55576];
      [1.81836,-0.95112,3.51122,2.92416,4.23544,4.43524];
};
% p-d
WS2_MoSe2_lay1_zpd={4.2090,4.3090,4.4090,4.5090,4.6090,4.7090,4.8090,4.9090,5.0090};
WS2_MoSe2_pd_lay1_inter={
      [-126.95384,-608231.81490,-4.78790,-9.54648,2.37464,1.73456,-4.81590,-2.55593];
      [-63.62313,-456752.02627,-4.32641,-9.30464,2.22383,1.72313,-4.20103,-2.61827];
      [-60.10128,-416354.59623,-4.25986,-9.18679,2.14489,1.70876,-3.95695,-2.66745];
      [-88.33519,-374410.49961,-4.47167,-9.06616,2.11711,1.69358,-3.98450,-2.70847];
      [-135.45587,-357616.99354,-4.70815,-8.98624,2.09580,1.67574,-4.04809,-2.73561];
      [-208.87368,-359445.05646,-4.94670,-8.94054,2.07682,1.65528,-4.12527,-2.74849];
      [-324.88401,-362265.13088,-5.18967,-8.90066,2.05890,1.62916,-4.20740,-2.72886];
      [-505.97833,-347068.43915,-5.43212,-8.83413,2.04359,1.60130,-4.30403,-2.69451];
      [-783.33619,-326848.72705,-5.66875,-8.76043,2.03028,1.57484,-4.41218,-2.66164];
};
WS2_MoSe2_lay2_zpd={4.2090,4.3090,4.4090,4.5090,4.6090,4.7090,4.8090,4.9090,5.0090};
WS2_MoSe2_pd_lay2_inter={
      [669.69429,60394.77167,-5.78887,-8.16467,2.15800,1.54469,-3.61924,-1.26630];
      [753.37748,47216.60979,-5.82092,-7.95757,2.08852,1.49695,-3.43250,-1.12206];
      [1102.20355,43971.09243,-6.01775,-7.86023,2.06840,1.46517,-3.49568,-1.06103];
      [1402.28949,44375.12177,-6.12587,-7.81734,2.02450,1.43461,-3.43800,-1.00617];
      [1920.36651,52390.28217,-6.28089,-7.87096,1.98922,1.41950,-3.42295,-1.04058];
      [2635.27316,63858.63115,-6.43633,-7.94514,1.95911,1.40639,-3.43270,-1.08635];
      [3674.01826,78321.91499,-6.60140,-8.02427,1.93246,1.39259,-3.45876,-1.12746];
      [5263.07092,101407.66323,-6.78315,-8.13662,1.90866,1.38386,-3.49768,-1.19967];
      [7115.95251,237376.84448,-6.93046,-8.59220,1.85675,1.43091,-3.37280,-1.62762];
};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WS2-WSe2
% p-p
WS2_WSe2_zpp={0.0,2.6180,2.7180,2.8180,2.9180,3.0180,3.1180,3.2180,3.3180};
WS2_WSe2_pp_inter={
      [0.00000, 0.00000,0.00000,0.00000,0.00000,0.00000];
      [3.92000,-1.57602,2.93276,2.57963,3.74489,5.38711];
      [3.48936,-1.22990,3.02359,2.70567,3.86727,5.74245];
      [3.15243,-1.33089,3.10122,2.69872,3.94262,5.29563];
      [2.81221,-1.39387,3.18922,2.70165,4.04606,4.95728];
      [2.51569,-1.43231,3.27344,2.70799,4.13632,4.68777];
      [2.25836,-1.52720,3.35505,2.69566,4.21617,4.39146];
      [2.00990,-1.67156,3.44402,2.66578,4.30962,4.06663];
      [1.79743,-2.49663,3.52780,2.48173,4.39094,3.44422];
};
% p-d
WS2_WSe2_lay1_zpd={0.0,4.2090,4.3090,4.4090,4.5090,4.6090,4.7090,4.8090,4.9090};
WS2_WSe2_pd_lay1_inter={
      [0.00000, 0.0000000000, 0.00000, 0.00000,0.00000,0.00000, 0.00000, 0.00000];
      [-295.56533,-552535.88580,-5.30125,-9.43315,2.51127,1.78879,-5.56547,-2.82712];
      [-93.38387,-391146.41886,-4.55082,-9.15099,2.35236,1.79349,-4.93861,-2.97764];
      [-95.96507,-220024.05370,-4.54377,-8.73929,2.30049,1.79307,-4.86140,-3.09680];
      [-107.33666,-170455.05744,-4.59007,-8.53352,2.26301,1.78298,-4.85663,-3.17982];
      [-136.12446,-176585.12863,-4.71576,-8.50672,2.23781,1.78035,-4.91954,-3.31021];
      [-170.85097,-127382.90478,-4.83512,-8.27151,2.21535,1.75925,-4.99581,-3.33175];
      [-221.43106,-131902.68379,-4.97169,-8.25339,2.19804,1.74735,-5.09896,-3.40844];
      [-298.38868,-110448.77637,-5.13264,-8.11476,2.18238,1.72522,-5.21187,-3.42119];
};
WS2_WSe2_lay2_zpd={0.0,4.2090,4.3090,4.4090,4.5090,4.6090,4.7090,4.8090,4.9090};
WS2_WSe2_pd_lay2_inter={
      [0.00000, 0.0000000000, 0.00000, 0.00000,0.00000,0.00000, 0.00000, 0.00000];
      [309.62957,19870.07036,-5.28607,-7.45918,2.11860,1.52113,-3.53575,-1.15286];
      [159.67270,2224.96764,-4.80987,-6.06828,1.98459,1.28461,-2.97925,0.06975];
      [188.76768,3058.89528,-4.87838,-6.21595,1.94887,1.26772,-2.94149,0.07679];
      [311.63195,3984.88932,-5.15680,-6.33284,1.95437,1.24635,-3.12421,0.10666];
      [386.59387,5119.99774,-5.25579,-6.44303,1.93291,1.22142,-3.15570,0.15179];
      [460.06310,8932.76304,-5.32532,-6.74391,1.91580,1.26369,-3.22328,-0.21084];
      [744.92056,10043.45440,-5.59040,-6.76859,1.90738,1.20621,-3.32606,0.03502];
      [1104.60104,14664.60882,-5.79788,-6.95713,1.88843,1.21953,-3.38671,-0.14843];
};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MoSe2-WSe2
% p-p
MoSe2_WSe2_zpp={2.4850,2.5850,2.6850,2.7850,2.8850,2.9850,3.0850,3.1850,3.2850};
MoSe2_WSe2_pp_inter={
      [4.65311,-1.79133,2.87302,2.55163,3.57415,5.30139];
      [4.33977,-1.46749,2.93041,2.66289,3.61438,5.64587];
      [3.98009,-1.26448,2.99818,2.74881,3.67592,5.79404];
      [3.64490,-1.11078,3.06592,2.82628,3.73220,5.86999];
      [3.35169,-0.97766,3.13048,2.90311,3.77453,5.91835];
      [3.08740,-0.88039,3.19271,2.96857,3.80643,5.86950];
      [2.82712,-0.82439,3.25957,3.01411,3.84737,5.67251];
      [2.58675,-0.84085,3.32638,3.01450,3.88484,5.23662];
      [2.38127,-1.50940,3.38847,2.71856,3.90453,3.86883];
};
% p-d
MoSe2_WSe2_lay1_zpd={4.1425,4.2425,4.3425,4.4425,4.5425,4.6425,4.7425,4.8425,4.9425};
MoSe2_WSe2_pd_lay1_inter={
      [-108011.89278,-3481700.96796,-9.10365,-10.63900,2.73143,1.52861,-6.49996,-1.60049];
      [-58.19924,-1775308.02237,-4.23450,-10.15138,2.06072,1.55652,-3.28777,-1.80973];
      [-71.20572,-1475988.84430,-4.33756,-9.96980,2.01101,1.56459,-3.18302,-1.94567];
      [-91.41953,-1110776.93471,-4.46713,-9.73308,1.98015,1.57344,-3.18340,-2.08577];
      [-134.65732,-894248.08273,-4.68294,-9.54377,1.96280,1.55860,-3.25429,-2.11463];
      [-198.39007,-700602.50238,-4.89698,-9.34311,1.94918,1.54094,-3.34739,-2.12728];
      [-300.09405,-659818.85416,-5.12459,-9.25774,1.93771,1.52771,-3.45399,-2.16358];
      [-444.52877,-565938.92707,-5.33746,-9.11989,1.93054,1.51615,-3.58963,-2.20338];
      [-707.45410,-496331.71688,-5.59274,-9.00256,1.91159,1.48140,-3.65332,-2.11647];
};
MoSe2_WSe2_lay2_zpd={4.1425,4.2425,4.3425,4.4425,4.5425,4.6425,4.7425,4.8425,4.9425};
MoSe2_WSe2_pd_lay2_inter={
      [25213.86809,335941.69891,-8.13706,-9.18276,2.46036,1.68045,-5.44913,-2.17660];
      [34259.79079,109688.13881,-8.27029,-8.42173,2.50225,1.62081,-5.78051,-1.98519];
      [1345.41360,145319.21824,-6.17945,-8.53365,2.28173,1.61394,-4.80334,-2.08919];
      [2149.09500,141400.36058,-6.43606,-8.46391,2.25250,1.59180,-4.81734,-2.08668];
      [2103.48119,105396.16009,-6.37915,-8.23283,2.19693,1.56510,-4.68299,-2.03513];
      [3306.61960,53058.72648,-6.62035,-7.77587,2.13434,1.52698,-4.52235,-1.90815];
      [3781.94973,34839.57519,-6.66263,-7.48838,2.07740,1.46837,-4.37099,-1.67561];
      [3278.31628,36292.72359,-6.53547,-7.47797,2.01184,1.43840,-4.16359,-1.62257];
      [4596.95852,55010.46314,-6.69469,-7.67647,1.96210,1.44979,-4.05251,-1.80459];
};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assign parameters as global cell arrays
if(multilayer)
    il = 0;
for ilayer = 1 : nlayers - 1
   % Interlayer parameters for each bilayer configuration
   if(tmdc(1+il)==1 && tmdc(2+il)==1)
       %p-p
       pp_vint_z(:,ilayer)=MoS2_MoS2_zpp;
       pp_vint_parm(:,ilayer)=MoS2_MoS2_pp_inter;
       % lay1 pd - bottom metal to top chalcogen
       pd_vint_lay1_z(:,ilayer)=MoS2_MoS2_lay1_zpd;
       pd_vint_lay1_parm(:,ilayer)=MoS2_MoS2_pd_lay1_inter;
       % lay2 pd - top metal to bottom chalcogen
       pd_vint_lay2_z(:,ilayer)=MoS2_MoS2_lay2_zpd;
       pd_vint_lay2_parm(:,ilayer)=MoS2_MoS2_pd_lay2_inter;
   elseif(tmdc(1+il)==1 && tmdc(2+il)==2)
       %pp
       pp_vint_z(:,ilayer)=MoS2_MoSe2_zpp;
       pp_vint_parm(:,ilayer)=MoS2_MoSe2_pp_inter;
       %pd
       pd_vint_lay1_z(:,ilayer)=MoS2_MoSe2_lay1_zpd;
       pd_vint_lay1_parm(:,ilayer)=MoS2_MoSe2_pd_lay1_inter;
       pd_vint_lay2_z(:,ilayer)=MoS2_MoSe2_lay2_zpd;
       pd_vint_lay2_parm(:,ilayer)=MoS2_MoSe2_pd_lay2_inter;        
   elseif(tmdc(1+il)==1 && tmdc(2+il)==3)
       %pp
       pp_vint_z(:,ilayer)=MoS2_WS2_zpp;
       pp_vint_parm(:,ilayer)=MoS2_WS2_pp_inter;
       %pd
       pd_vint_lay1_z(:,ilayer)=MoS2_WS2_lay1_zpd;
       pd_vint_lay1_parm(:,ilayer)=MoS2_WS2_pd_lay1_inter;
       pd_vint_lay2_z(:,ilayer)=MoS2_WS2_lay2_zpd;
       pd_vint_lay2_parm(:,ilayer)=MoS2_WS2_pd_lay2_inter;         
   elseif(tmdc(1+il)==1 && tmdc(2+il)==4)
       %pp
       pp_vint_z(:,ilayer)=MoS2_WSe2_zpp;
       pp_vint_parm(:,ilayer)=MoS2_WSe2_pp_inter;
       %pd
       pd_vint_lay1_z(:,ilayer)=MoS2_WSe2_lay1_zpd;
       pd_vint_lay1_parm(:,ilayer)=MoS2_WSe2_pd_lay1_inter;
       pd_vint_lay2_z(:,ilayer)=MoS2_WSe2_lay2_zpd;
       pd_vint_lay2_parm(:,ilayer)=MoS2_WSe2_pd_lay2_inter;         
   elseif(tmdc(1+il)==2 && tmdc(2+il)==1)
       %pp
       pp_vint_z(:,ilayer)=MoS2_MoSe2_zpp;
       pp_vint_parm(:,ilayer)=MoS2_MoSe2_pp_inter;
       %pd
       pd_vint_lay1_z(:,ilayer)=MoS2_MoSe2_lay2_zpd;
       pd_vint_lay1_parm(:,ilayer)=MoS2_MoSe2_pd_lay2_inter;
       pd_vint_lay2_z(:,ilayer)=MoS2_MoSe2_lay1_zpd;
       pd_vint_lay2_parm(:,ilayer)=MoS2_MoSe2_pd_lay1_inter;          
   elseif(tmdc(1+il)==2 && tmdc(2+il)==2)
       %pp
       pp_vint_z(:,ilayer)=MoSe2_MoSe2_zpp;
       pp_vint_parm(:,ilayer)=MoSe2_MoSe2_pp_inter;
       %pd
       pd_vint_lay1_z(:,ilayer)=MoSe2_MoSe2_lay1_zpd;
       pd_vint_lay1_parm(:,ilayer)=MoSe2_MoSe2_pd_lay1_inter;
       pd_vint_lay2_z(:,ilayer)=MoSe2_MoSe2_lay2_zpd;
       pd_vint_lay2_parm(:,ilayer)=MoSe2_MoSe2_pd_lay2_inter;          
   elseif(tmdc(1+il)==2 && tmdc(2+il)==3)
       %pp
       pp_vint_z(:,ilayer)=WS2_MoSe2_zpp;
       pp_vint_parm(:,ilayer)=WS2_MoSe2_pp_inter;
       %pd
       pd_vint_lay1_z(:,ilayer)=WS2_MoSe2_lay2_zpd;
       pd_vint_lay1_parm(:,ilayer)=WS2_MoSe2_pd_lay2_inter;
       pd_vint_lay2_z(:,ilayer)=WS2_MoSe2_lay1_zpd;
       pd_vint_lay2_parm(:,ilayer)=WS2_MoSe2_pd_lay1_inter;           
   elseif(tmdc(1+il)==2 && tmdc(2+il)==4)
       %pp
       pp_vint_z(:,ilayer)=MoSe2_WSe2_zpp;
       pp_vint_parm(:,ilayer)=MoSe2_WSe2_pp_inter;
       %pd
       pd_vint_lay1_z(:,ilayer)=MoSe2_WSe2_lay1_zpd;
       pd_vint_lay1_parm(:,ilayer)=MoSe2_WSe2_pd_lay1_inter;
       pd_vint_lay2_z(:,ilayer)=MoSe2_WSe2_lay2_zpd;
       pd_vint_lay2_parm(:,ilayer)=MoSe2_WSe2_pd_lay2_inter;         
   elseif(tmdc(1+il)==3 && tmdc(2+il)==1)
       %pp
       pp_vint_z(:,ilayer)=MoS2_WS2_zpp;
       pp_vint_parm(:,ilayer)=MoS2_WS2_pp_inter;
       %pd
       pd_vint_lay1_z(:,ilayer)=MoS2_WS2_lay2_zpd;
       pd_vint_lay1_parm(:,ilayer)=MoS2_WS2_pd_lay2_inter;
       pd_vint_lay2_z(:,ilayer)=MoS2_WS2_lay1_zpd;
       pd_vint_lay2_parm(:,ilayer)=MoS2_WS2_pd_lay1_inter;       
   elseif(tmdc(1+il)==3 && tmdc(2+il)==2)
       %pp
       pp_vint_z(:,ilayer)=WS2_MoSe2_zpp;
       pp_vint_parm(:,ilayer)=WS2_MoSe2_pp_inter;
       %pd
       pd_vint_lay1_z(:,ilayer)=WS2_MoSe2_lay1_zpd;
       pd_vint_lay1_parm(:,ilayer)=WS2_MoSe2_pd_lay1_inter;
       pd_vint_lay2_z(:,ilayer)=WS2_MoSe2_lay2_zpd;
       pd_vint_lay2_parm(:,ilayer)=WS2_MoSe2_pd_lay2_inter;     
   elseif(tmdc(1+il)==3 && tmdc(2+il)==3)
       %pp
       pp_vint_z(:,ilayer)=WS2_WS2_zpp;
       pp_vint_parm(:,ilayer)=WS2_WS2_pp_inter;
       %pd
       pd_vint_lay1_z(:,ilayer)=WS2_WS2_lay1_zpd;
       pd_vint_lay1_parm(:,ilayer)=WS2_WS2_pd_lay1_inter;
       pd_vint_lay2_z(:,ilayer)=WS2_WS2_lay2_zpd;
       pd_vint_lay2_parm(:,ilayer)=WS2_WS2_pd_lay2_inter; 
   elseif(tmdc(1+il)==3 && tmdc(2+il)==4)
       %pp
       pp_vint_z(:,ilayer)=WS2_WSe2_zpp;
       pp_vint_parm(:,ilayer)=WS2_WSe2_pp_inter;
       %pd
       pd_vint_lay1_z(:,ilayer)=WS2_WSe2_lay1_zpd;
       pd_vint_lay1_parm(:,ilayer)=WS2_WSe2_pd_lay1_inter;
       pd_vint_lay2_z(:,ilayer)=WS2_WSe2_lay2_zpd;
       pd_vint_lay2_parm(:,ilayer)=WS2_WSe2_pd_lay2_inter;        
   elseif(tmdc(1+il)==4 && tmdc(2+il)==1)
       %pp
       pp_vint_z(:,ilayer)=MoS2_WSe2_zpp;
       pp_vint_parm(:,ilayer)=MoS2_WSe2_pp_inter;
       %pd
       pd_vint_lay1_z(:,ilayer)=MoS2_WSe2_lay2_zpd;
       pd_vint_lay1_parm(:,ilayer)=MoS2_WSe2_pd_lay2_inter;
       pd_vint_lay2_z(:,ilayer)=MoS2_WSe2_lay1_zpd;
       pd_vint_lay2_parm(:,ilayer)=MoS2_WSe2_pd_lay1_inter;       
   elseif(tmdc(1+il)==4 && tmdc(2+il)==2)
       %pp
       pp_vint_z(:,ilayer)=MoSe2_WSe2_zpp;
       pp_vint_parm(:,ilayer)=MoSe2_WSe2_pp_inter;
       %pd
       pd_vint_lay1_z(:,ilayer)=MoSe2_WSe2_lay2_zpd;
       pd_vint_lay1_parm(:,ilayer)=MoSe2_WSe2_pd_lay2_inter;
       pd_vint_lay2_z(:,ilayer)=MoSe2_WSe2_lay1_zpd;
       pd_vint_lay2_parm(:,ilayer)=MoSe2_WSe2_pd_lay1_inter;  
   elseif(tmdc(1+il)==4 && tmdc(2+il)==3)
       %pp
       pp_vint_z(:,ilayer)=WS2_WSe2_zpp;
       pp_vint_parm(:,ilayer)=WS2_WSe2_pp_inter;
       %pd
       pd_vint_lay1_z(:,ilayer)=WS2_WSe2_lay2_zpd;
       pd_vint_lay1_parm(:,ilayer)=WS2_WSe2_pd_lay2_inter;
       pd_vint_lay2_z(:,ilayer)=WS2_WSe2_lay1_zpd;
       pd_vint_lay2_parm(:,ilayer)=WS2_WSe2_pd_lay1_inter; 
   elseif(tmdc(1+il)==4 && tmdc(2+il)==4)     
       %pp
       pp_vint_z(:,ilayer)=WSe2_WSe2_zpp;
       pp_vint_parm(:,ilayer)=WSe2_WSe2_pp_inter;
       %pd
       pd_vint_lay1_z(:,ilayer)=WSe2_WSe2_lay1_zpd;
       pd_vint_lay1_parm(:,ilayer)=WSe2_WSe2_pd_lay1_inter;
       pd_vint_lay2_z(:,ilayer)=WSe2_WSe2_lay2_zpd;
       pd_vint_lay2_parm(:,ilayer)=WSe2_WSe2_pd_lay2_inter;        
   end
   il = il + 1;
end
else
    % Not multilayer - so set empty parameters
    pp_vint_z={};
    pp_vint_parm={};
    pd_vint_lay1_z={};
    pd_vint_lay1_parm={};
    pd_vint_lay2_z={};
    pd_vint_lay2_parm={};
end
end

