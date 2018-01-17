%:::::::::::::::::::::% 
%                     %
% Physical properties %
%                     %
%:::::::::::::::::::::%
  1
  1   FLUID  .01  1.00  

%:::::::::::::::::::::% 
%                     %
% Boundary conditions %
%                     %
%:::::::::::::::::::::%
   4 
%---------------------------------------------------%
% This is the default boundary condition, allways 1 %  
%---------------------------------------------------%
   1         WALL     0.0      0.0   0.0
%========%
% INFLOW %
%========%
   2         INFLOW    FILE  parabolic-0-2.prof 
%=========%
% OUTFLOW %
%=========%
   3         OUTFLOW   FILE  parabolic-0-2.prof  
%==========%
% SYMMETRY %
%==========%
   4         SYMMETRY  0.0      0.0    0.0 

%::::::::::::::::::::% 
%                    %
% Initial conditions %
%                    %
%::::::::::::::::::::%
   1
   1   1.0   0.0   0.0
   
