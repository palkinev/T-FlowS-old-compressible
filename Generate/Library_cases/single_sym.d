#======================================================================#
#
#                               GO WITH THE FLOW !
#                       
#                45------46------47------------------------------48
#                /|                                              /|
#              41 |                                             / |
#              /  |                                            /  |
#            37   |                                           /   |
#            /   29                                          /   32
#   Flow   33------34------35------------------------------36     |
#  =---->   |     |    26------27                           |     |
#           |    13----/-14----/|15-------------------------|----16  
#           |    /   22------23 |                           |    /
#           |   9     |10-----|11                           |  12 
#          17  /      |/      |/                           20  /
#           | 5       6-------7                             | 8
#           |/                                              |/
#           1-------2-------3-------------------------------4
#                             
#----------------------------------------------------------------------#

#===========================================#
#  Nodes (cells), boundary cells and sides  #
#===========================================#
  30000 10000 90000

#========#
# Points #
#========#
48
#--------------#
# ground floor #
#--------------#
 1  0.0   0.0    0.0
 2  1.00  0.0    0.0
 3  2.00  0.0    0.0
 4  6.00  0.0    0.0
 5  0.0   1.00   0.0
 6  1.00  1.00   0.0
 7  2.00  1.00   0.0
 8  6.00  1.00   0.0
 9  0.0   2.00   0.0
10  1.00  2.00   0.0
11  2.00  2.00   0.0
12  6.00  2.00   0.0
13  0.0   3.00   0.0
14  1.00  3.00   0.0
15  2.00  3.00   0.0
16  6.00  3.00   0.0
#--------------#
# middle floor #
#--------------#
17  0.0   0.0   1.000
18  1.00  0.0   1.000
19  2.00  0.0   1.000
20  6.00  0.0   1.000
21  0.0   1.00  1.000
22  1.00  1.00  1.000
23  2.00  1.00  1.000
24  6.00  1.00  1.000
25  0.0   2.00  1.000
26  1.00  2.00  1.000
27  2.00  2.00  1.000
28  6.00  2.00  1.000
29  0.0   3.00  1.000
30  1.00  3.00  1.000
31  2.00  3.00  1.000
32  6.00  3.00  1.000
#-----------#
# top floor #
#-----------#
33  0.0   0.0   2.00
34  1.00  0.0   2.00
35  2.00  0.0   2.00
36  6.00  0.0   2.00
37  0.0   1.00  2.00
38  1.00  1.00  2.00
39  2.00  1.00  2.00
40  6.00  1.00  2.00
41  0.0   2.00  2.00
42  1.00  2.00  2.00
43  2.00  2.00  2.00
44  6.00  2.00  2.00
45  0.0   3.00  2.00
46  1.00  3.00  2.00
47  2.00  3.00  2.00
48  6.00  3.00  2.00
#========#
# Blocks #
#========#
17
#---------------------------------------
#                 13---14---15---16
#                  |  3 |  5 |  7 |
#                  9---10---11---12
# first floor      |  2 |    |  8 |
#                  5----6----7----8
#                  |  1 |  4 |  6 |
#                  1----2----3----4
#---------------------------------------
   1   8  8  8
       1.0  1.0  1.0
       1  2  5  6 17 18 21 22
   2   8  8  8
       1.0  1.0  1.0
       5  6  9 10 21 22 25 26
   3   8  8  8
       1.0  1.0  1.0
       9 10 13 14 25 26 29 30
   4    8  8  8
       1.0  1.0  1.0
       2  3  6  7 18 19 22 23
   5    8  8  8
       1.0  1.0  1.0
      10 11 14 15 26 27 30 31
   6   15  8  8
       .25  1.0  1.0
       3  4  7  8 19 20 23 24
   7   15  8  8
       .25  1.0  1.0
      11 12 15 16 27 28 31 32
   8   15  8  8
       .25  1.0  1.0
       7  8 11 12 23 24 27 28
#---------------------------------------
#                 29---30---31---32
#                  | 11 | 13 | 15 |
#                 25---26---27---28
# second floor     | 10 | 17 | 16 |
#                 21---22---23---24
#                  |  9 | 12 | 14 |
#                 17---18---19---20
#---------------------------------------
   9   8  8  8
       1.0  1.0  1.0
      17 18 21 22 33 34 37 38
  10   8  8  8
       1.0  1.0  1.0
      21 22 25 26 37 38 41 42
  11   8  8  8
       1.0  1.0  1.0
      25 26 29 30 41 42 45 46
  12    8  8  8
       1.0  1.0  1.0
      18 19 22 23 34 35 38 39
  13    8  8  8
       1.0  1.0  1.0
      26 27 30 31 42 43 46 47
  14  15  8  8 
       .25  1.0  1.0
      19 20 23 24 35 36 39 40
  15  15  8  8
       .25  1.0  1.0
      27 28 31 32 43 44 47 48
  16  15  8  8
       .25  1.0  1.0
      23 24 27 28 39 40 43 44
  17    8  8  8
       1.0  1.0  1.0
      22 23 26 27 38 39 42 43
%%%%%%%%%%%%%%%%%%%%
%   line section   %
%%%%%%%%%%%%%%%%%%%%
  0
%%%%%%%%%%%%%%%%%%%%%%%
%   surface section   %  -> still under construction
%%%%%%%%%%%%%%%%%%%%%%%
  0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   boundary conditions and material data   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 24 
#--------#
# Inflow #
#--------#
    1     1  1  1     1  7  7
       1  2
    2     1  1  1     1  7  7
       2  2
    3     1  1  1     1  7  7
       3  2
    4     1  1  1     1  7  7
       9  2
    5     1  1  1     1  7  7
      10  2
    6     1  1  1     1  7  7
      11  2
#---------#
# Outflow #
#---------#
    7    14  1  1    14  7  7
       6  3
    8    14  1  1    14  7  7
       7  3
    9    14  1  1    14  7  7
       8  3
   10    14  1  1    14  7  7
      14  3
   11    14  1  1    14  7  7
      15  3
   12    14  1  1    14  7  7
      16  3
#----------#
# Symmetry #
#----------#
% lower blocks 
   13     1  1  1     7  1  7
       1  4
   14     1  7  1     7  7  7
       3  4
   15     1  1  1     7  1  7
       4  4
   16     1  7  1     7  7  7
       5  4
   17     1  1  1    14  1  7
       6  4
   18     1  7  1    14  7  7
       7  4
% upper blocks 
   19     1  1  1     7  1  7
       9  4
   20     1  7  1     7  7  7
      11  4
   21     1  1  1     7  1  7
      12  4
   22     1  7  1     7  7  7
      13  4
   23     1  1  1    14  1  7
      14  4
   24     1  7  1    14  7  7
      15  4
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   periodic boundaries   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
  0
%%%%%%%%%%%%%%%%%%
%   refinement   %
%%%%%%%%%%%%%%%%%%
  0
%%%%%%%%%%%%%%%%%%%%
%   don't smooth   %
%%%%%%%%%%%%%%%%%%%%
  1
     1
     20  0.9   
     0.0  0.0  0.0  6.0  3.0  2.0












