%/////////////////////%
% Phisical properties %
%/////////////////////%
   2
   1  FLUID .005  1.00  
   2  SOLID .005  1.00  

%/////////////////////%
% Boundary conditions %
%/////////////////////%
   3 
%---------------------------------------------------%
% This is the default boundary condition, allways 1 %  
%---------------------------------------------------%
   1        WALL       0.0      0.0   0.0
%========%
% INFLOW %
%========%
   2        INFLOW    FILE  parabolic-0-2.prof 
%=========%
% OUTFLOW %
%=========%
   3        OUTFLOW   FILE  parabolic-0-2.prof 
   
%////////////////////%
% Initial conditions %
%////////////////////%
   2
   1   1.0   0.0   0.0   0.0
   2   0.0   0.0   0.0   1.0
