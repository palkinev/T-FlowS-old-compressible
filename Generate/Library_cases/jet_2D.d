#
#
#                15----------------------16
#                /|                      /|
#              13----------------------14 |
#               | |         (3)         | |
#    10---------|11---------------------|12
#    /|         |/|                     |/|
#   7-----------8-----------------------9 |
#   | |         | |                     | |
#   | |         | |                     | |
#   | |   (1)   | |         (2)         | |
#   | |         | |                     | |
#   | |         | |                     | |
#   | 4 - - - - | 5- - - - - - - - - - -|-6
#   |/          |/                      |/ 
#   1-----------2-----------------------3  
#
#-------------------------------------------#
#  Nodes (cells), boundary cells and sides  #
#-------------------------------------------#
  40000 13333 120000

#----------#
#  Points  #
#----------#
16
  1  0.0 0.0 0.0
  2  1.0 0.0 0.0
  3  9.0 0.0 0.0
  4  0.0 1.0 0.0
  5  1.0 1.0 0.0
  6  9.0 1.0 0.0
  7  0.0 0.0 2.0
  8  1.0 0.0 2.0
  9  9.0 0.0 2.0
 10  0.0 1.0 2.0
 11  1.0 1.0 2.0
 12  9.0 1.0 2.0
 13  1.0 0.0 3.0
 14  9.0 0.0 3.0
 15  1.0 1.0 3.0
 16  9.0 1.0 3.0
#----------#
#  Blocks  #
#----------#
3
  1   41   4  41 
      1.0  1.0  1.0 
      1  2  4  5  7  8 10 11
  2  121   4  41 
      0.5  1.0  1.0 
      2  3  5  6  8  9 11 12
  3  121   4  21 
      0.5  1.0  1.0 
      8  9 11 12 13 14 15 16
#--------#
#        #  
#--------#
   0 
   0
#-----------------------#
#  Boundary conditions  #
# (it will use default) #
#-----------------------#
  5
    1     Kmax
        1   2
    2     Imin
        1   3
#-----------------
    3     Imax
        2   4
#-----------------
    4     Imax
        3   4
#-----------------
    5     Kmax
        3   5
#-----------------------#
#  Periodic boundaries  #
#-----------------------#
  3
      1    1  2  8  7
           4  5 11 10
      2    2  3  9  8
           5  6 12 11
      3    8  9 14 13
          11 12 16 15
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