%/////////////////////%
% Phisical properties %
%/////////////////////%
   2
   1  SOLID .05  1.00  .005
   2  FLUID .05  1.00  .5 

%/////////////////////%
% Boundary conditions %
%/////////////////////%
   3 
%-- This is the default boundary condition, allways 1 %  
   1        Symmetry   0.0      0.0   0.0   0.0
%-- Left
   2        Wall       0.0      0.0   0.0   2.0
%-- Right 
   3        Wall       0.0      0.0   0.0   1.0 
   
%////////////////////%
% Initial conditions %
%////////////////////%
   2
   1   0.0   0.0   0.0   1.0
   2   0.0   0.0   0.0   2.0
