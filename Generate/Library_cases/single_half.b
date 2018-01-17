%:::::::::::::::::::::% 
%                     %
% PHYSICAL PROPERTIES %
%                     %
%:::::::::::::::::::::%
  1  
  1   fluid   .01  1.00  

%:::::::::::::::::::::% 
%                     %
% Boundary conditions %
%                     %
%:::::::::::::::::::::%
   4 
%---------------------------------------------------%
% This is the default boundary condition, allways 1 %  
%---------------------------------------------------%
   1         wall      0.0   0.0   0.0
%========%
% INFLOW %
%========%
   2         inflow    file  parabolic-0-2.prof 
%=========%
% OUTFLOW %
%=========%
   3         outflow   file  parabolic-0-2.prof 
%==========%
% SYMMETRY %
%==========%
   4         symmetry  0.0   0.0   0.0 
   
%::::::::::::::::::::% 
%                    %
% Initial conditions %
%                    %
%::::::::::::::::::::%
   1
   1   1.0   0.0   0.0
