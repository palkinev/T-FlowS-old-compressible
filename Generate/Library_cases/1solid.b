%/////////////////////%
% Phisical properties %
%/////////////////////%
   1
   1  SOLID .05  1.00  .1

%/////////////////////%
% Boundary conditions %
%/////////////////////%
   3 
%-- This is the default boundary condition, allways 1 %  
   1        Symmetry   0.0      0.0   0.0   0.0
%-- Left
   2        Wall       0.0      0.0   0.0   1.0
%-- Right 
   3        Wall       0.0      0.0   0.0   1.5 
   
%////////////////////%
% Initial conditions %
%////////////////////%
   1
   1   0.0   0.0   0.0   1.0
