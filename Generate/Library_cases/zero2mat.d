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
  1  -0.1             0.0            -0.1
  2   4.1             0.0            -0.1
  3  -0.1             0.1            -0.1
  4   4.1             0.1            -0.1
  5  -0.1             0.0             1.1            
  6   4.1             0.0             1.1            
  7  -0.1             0.1             1.1            
  8   4.1             0.1             1.1            
#----------#
#  Blocks  #
#----------#
1
  1  127  3  37 
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
  1
     1  4  1  4   123  2  33 
        1  2
#-----------------------#
#  Periodic boundaries  #
#-----------------------#
  1
      1    1 2 6 5
           3 4 8 7
#-------------------#
#  Copy boundaries  #
#-------------------#
   0
#------------
# Refinement
#-----------
   0
#------------
# Smoothing 
#-----------
   0

