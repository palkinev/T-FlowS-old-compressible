%-----------------------%
%                       %
%  PHYSICAL PROPERTIES  %
%                       %
%-----------------------%
  2
  1   FLUID   .005  1.00  
  2   FLUID   .005  1.00  

%-----------------------%
%                       %
%  BOUNDARY CONDITIONS  %
%                       %
%-----------------------%
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
   
%----------------------%
%                      %
%  INITIAL CONDITIONS  %
%                      %
%----------------------%
  2
  1   1.0  0.0  0.0
  2   1.0  0.0  0.0
