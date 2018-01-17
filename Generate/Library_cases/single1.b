%/////////////////////%
% Physical properties %
%/////////////////////%
  1
  1   FLUID  .005  1.00  
%/////////////////////%
% Boundary conditions %
%/////////////////////%
   6 
%---------------------------------------------------%
% This is the default boundary condition, allways 1 %  
%---------------------------------------------------%
   1        WALL       0.0      0.0   0.0
   2        WALL       0.0      0.0   0.0
   3        WALL       0.0      0.0   0.0
%========%
% INFLOW %
%========%
   4        INFLOW    FILE  parabolic-0-2.prof 
%=========%
% OUTFLOW %
%=========%
   5        OUTFLOW   FILE  parabolic-0-2.prof 
%=========%
% OUTFLOW %
%=========%
   6        SYMMETRY  FILE  parabolic-0-2.prof 
%////////////////////%
% Initial conditions %
%////////////////////%
  1
  1   1.0  0.0  0.0
   
