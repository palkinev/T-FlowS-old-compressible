%-------------------------------------------------------------%
%  Physical properties                                        % 
%-------------------------------------------------------------%
  4
       1  FLUID  0.01   1.00    0.02   1.0
       2  FLUID  0.01   1.00    0.02   1.0
       3  FLUID  0.01   1.00    0.02   1.0
       4  SOLID  0.01   1000.0  2.0    1.0
%-------------------------------------------------------------%
%  Boundary conditions                                        % 
%-------------------------------------------------------------%
   5 
       1    wall    0.0      0.0   0.0   0.0
       2    inflow  0.25     0.0   0.0   0.0
       3    inflow  0.0      0.0  -0.5   0.0
       4    outflow 0.5      0.0   0.0   0.0
       5    wall    0.0      0.0   0.0   1.0
%-------------------------------------------------------------%
%  Initial conditions                                         % 
%-------------------------------------------------------------%
   4 
       1   0.25 0.0  0.0   0.0
       2   0.0  0.0 -0.5   0.0
       3   0.25 0.0  0.0   0.0
       4   0.0  0.0  0.0   0.0
