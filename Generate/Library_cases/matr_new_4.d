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
   1   33  33  33 
      5.0  1.0  1.0
       1  9  3 11  5 13  7 15
   2   33  33  33 
      1.0  5.0  1.0
       1  2  9 10  5  6 13 14
   3   33  33  33   
      1.0  0.2  1.0
      11 12  3  4 15 16  7  8
   4   33  33  33   
      0.2  1.0  1.0
      10  2 12  4 14  6 16  8

  5   33  33  33   
      1.0  1.0  0.2 
      13 14 15 16  5  6  7  8

#---- Block on the ceiling
  6   33  33  17   
      1.0  1.0  5.0 
       5  6  7  8 17 18 19 20 

#---- Blocks on the floor
   7   33  33  17
      5.0  1.0  1.0 
      21 25 23 27  1  9  3 11
   8   33  33  17
      1.0  5.0  1.0 
      21 22 25 26  1  2  9 10
   9   33  33  17
      1.0  0.2  1.0 
      27 28 23 24 11 12  3  4
  10   33  33  17
      0.2  1.0  1.0 
      26 22 28 24 10  2 12  4
#==================#
# Lines & surfaces #
#==================#
  0
  0
#---------------------#
# Boundary conditions #
#---------------------#
   0
#-----------------------#
#  Periodic boundaries  #
#-----------------------#
   6
     1    1  3  7  5
          2  4  8  6
     2    1  2  6  5
          3  4  8  7
     3    5  7 19 17
          6  8 20 18
     4    5  6 18 17
          7  8 20 19
     5   21 23  3  1
         22 24  4  2
     6   21 22  2  1
         23 24  4  3
#------------#
# Refinement #
#------------#
   0
#-----------#
# Smoothing #
#-----------#
   2
     1  x y z  
       10  0.9
       0.0  0.0  0.0  1.2  1.2  1.0 
     2 
       1  0.9
       0.45  0.45  0.0  0.75  0.75  0.3 


