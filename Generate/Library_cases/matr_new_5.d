#===========================================
#
#         19-------------------20
#         /|                   /|
#        / |                  / |
#       /  7-----------------/--8
#      /  /|      6         /  /|
#     /  / |               /  / |
#   17-------------------18  /  |
#    | /   |      5       | /   |
#    |/    |              |/    |
#    5--------------------6     |
#    |     |    15---16   |     |
#    |     3--- /|   /| --|-----4
#    |    /   13---14 |   |    /|
#    |   /     |11--|12   |   / | 
#    |  / 23   |    |/    |4 / 24
#    | /       9---10     | /  /     
#    |/  /      27 - 28   |/  /10
#    1--------------------2  /
#    |    7   25---26     | /
#    |/         8         |/
#   21-------------------22
#                
#
#===========================================

#===========================================#
#  Nodes (cells), boundary cells and sides  #
#===========================================#
  30000 10000 90000

#========#
# Points #
#========#
28

/* domain */
 1  0.0  0.0  0.1
 2  1.2  0.0  0.1
 3  0.0  1.2  0.1
 4  1.2  1.2  0.1
 5  0.0  0.0  0.8
 6  1.2  0.0  0.8
 7  0.0  1.2  0.8
 8  1.2  1.2  0.8

/* cube */
 9  0.45 0.45 0.1
10  0.75 0.45 0.1
11  0.45 0.75 0.1
12  0.75 0.75 0.1
13  0.45 0.45 0.3
14  0.75 0.45 0.3
15  0.45 0.75 0.3
16  0.75 0.75 0.3

/* top */
17  0.0  0.0  1.0
18  1.2  0.0  1.0
19  0.0  1.2  1.0
20  1.2  1.2  1.0

/* floor */
21  0.0  0.0  0.0
22  1.2  0.0  0.0
23  0.0  1.2  0.0
24  1.2  1.2  0.0
25  0.45 0.45 0.0
26  0.75 0.45 0.0
27  0.45 0.75 0.0
28  0.75 0.75 0.0

#========#
# Blocks #
#========#
10
#---- Blocks around the cube
   1    9   9   9 
      5.0  1.0  1.0
       1  9  3 11  5 13  7 15
   2    9   9   9 
      1.0  5.0  1.0
       1  2  9 10  5  6 13 14
   3    9   9   9   
      1.0  0.2  1.0
      11 12  3  4 15 16  7  8
   4    9   9   9   
      0.2  1.0  1.0
      10  2 12  4 14  6 16  8

#---- Block on the top of the cube
  5    9   9   9   
      1.0  1.0  0.2 
      13 14 15 16  5  6  7  8

#---- Block on the ceiling
  6    9   9   5   
      1.0  1.0  5.0 
       5  6  7  8 17 18 19 20 

#---- Blocks on the floor
   7    9   9   5
      5.0  1.0  0.2 
      21 25 23 27  1  9  3 11
   8    9   9   5
      1.0  5.0  0.2 
      21 22 25 26  1  2  9 10
   9    9   9   5
      1.0  0.2  0.2 
      27 28 23 24 11 12  3  4
  10    9   9   5
      0.2  1.0  0.2 
      26 22 28 24 10  2 12  4

#=======#
# Lines # 
#=======#
 16 
# on the top of the cube
    -1  13  14
        -0.92
    -2  15  16
        -0.92
    -3  13  15
        -0.92
    -4  14  16
        -0.92
# ..................
    -5   9  13
        5.0
    -6  10  14
        5.0
    -7  11  15
        5.0
    -8  12  16
        5.0
# lines between upper and lower cube surfaces  
    -9   9  10
        -0.92
   -10  11  12 
        -0.92
   -11   9  11
        -0.92
   -12  10  12 
        -0.92
# lines between upper and lower cube surfaces  
   -13  25  26
        -0.92
   -14  27  28 
        -0.92
   -15  25  27
        -0.92
   -16  26  28 
        -0.92
#==========#
# Surfaces #
#==========#
  5
# surface on the top of the cube
     1  13 14 16 15
        -0.92 1. 1.
# surfaces on the upper sides of the cubes
     2   9 11 15 14
        1.0 1.0 5.0
     3  10 12 16 14
        1.0 1.0 5.0
     4   9 10 14 13
        1.0 1.0 5.0
     5  11 12 16 15
        1.0 1.0 5.0
#---------------------#
# Boundary conditions #
#---------------------#
   0
#-----------------------#
#  Periodic boundaries  #
#-----------------------#
   0
#   6
#     1    1  3  7  5
#          2  4  8  6
#     2    1  2  6  5
#          3  4  8  7
#     3    5  7 19 17
#          6  8 20 18
#     4    5  6 18 17
#          7  8 20 19
#     5   21 23  3  1
#         22 24  4  2
#     6   21 22  2  1
#         23 24  4  3
#------------#
# Refinement #
#------------#
   0
#-----------#
# Smoothing #
#-----------#
   1
     1  x y z  
       0  0.
       0.0  0.0  0.0  1.2  1.2  1.0 



