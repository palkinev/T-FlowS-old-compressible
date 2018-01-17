#
#
#     7-----------8
#    /|          /|
#   5-----------6 |
#   | |         | |
#   | |         | |
#   | |         | |
#   | 3 - - - - | 4
#   |/          |/
#   1-----------2
#
#

#-------------------------------------------#
#  Nodes (cells), boundary cells and sides  #
#-------------------------------------------#
  30000 10000 90000

#----------#
#  Points  #
#----------#
8
  1   0.0             0.0             0.0
  2   4.5             0.0             0.0
  3   0.0             1.0             0.0
  4   4.5             1.0             0.0
  5   0.0             0.0             1.0            
  6   4.5             0.0             1.0            
  7   0.0             1.0             1.0            
  8   4.5             1.0             1.0            
#----------#
#  Blocks  #
#----------#
1
  1   19  5  5 
      1.0  1.0  1.0 
      1  2  3  4  5  6  7  8
#--------#
#        #  
#--------#
   0 
   0
#-----------------------#
#  Boundary conditions  #
# (it will use default) #
#-----------------------#
   0
#-----------------------#
#  Periodic boundaries  #
#-----------------------#
   0
#------------
# Refinement
#-----------
   1
     1  1
        1  4.0  0.0  0.0  4.5  1.0  1.0 
#------------
# Smoothing 
#-----------
   1
     1
       0  0.0    
       0.0  0.0  0.0  1.0  1.0  1.0

