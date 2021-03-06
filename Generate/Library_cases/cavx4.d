%======================================================
%
%       16------------------17------------------18
%       /|                  /                   /| 
%     13------------------14------------------15 |
%      | |                                     | |     
%      |10 - - - - - - - - -11 - - - - - - - - |12     
%      |/|                  /                  |/|     
%      7 - - - - - - - - - 8 - - - - - - - - - 9 |     
%      | |                                     | |     
%      | 4-------------------5-----------------|-6  
%      |/                   /                  |/ 
%      1-------------------2-------------------3
%
%======================================================

#-------------------------------------------#
#  Nodes (cells), boundary cells and sides  #
#-------------------------------------------#
  30000 10000 90000

#----------#
#  Points  #
#----------#
18
# floor
  1   0.0  0.0  0.0
  2   0.1  0.0  0.0 
  3   1.0  0.0  0.0
  4   0.0  0.5  0.0
  5   0.1  0.5  0.0 
  6   1.0  0.5  0.0
# middle
  7   0.0  0.0  0.9 
  8   0.5  0.0  0.5 
  9   1.0  0.0  0.1
 10   0.0  0.5  0.9
 11   0.5  0.5  0.5 
 12   1.0  0.5  0.1
# roof 
 13   0.0  0.0  1.0
 14   0.9  0.0  1.0 
 15   1.0  0.0  1.0
 16   0.0  0.5  1.0
 17   0.9  0.5  1.0 
 18   1.0  0.5  1.0
#----------#
#  Blocks  #
#----------#
4
  1   11  4  11 
      1.0 1.0 1.0
#     -0.9  1.0  -0.9 
      1  2  4  5  7  8 10 11
  2   11  4  11 
      1.0 1.0 1.0
#     -0.9  1.0  -0.9 
      2  3  5  6  8  9 11 12
  3   11  4  11 
      1.0 1.0 1.0
#      -0.9  1.0  -0.9
      7  8 10 11 13 14 16 17
  4   11  4  11 
      1.0 1.0 1.0
#      -0.9  1.0  -0.9
      8  9 11 12 14 15 17 18
#---------#
#  Lines  #  
#---------#
   0
#--------#
#        #  
#--------#
   0
#--------#
#        #  
#--------#
   4
   1     1   1   10   10  3  10
         3   2
   2     1   1   10   10  3  10
         4   2
   3     1   1   1    10  3   1
         1   3
   4     1   1   1    10  3   1
         2   3
#-----------------------#
#  Periodic boundaries  #
#-----------------------#
   4  
   1   1  2  8  7
       4  5 11 10   
   2   2  3  9  8
       5  6 12 11   
   3   7  8 14 13
      10 11 17 16   
   4   8  9 15 14
      11 12 18 17   
#--------------#
#  Refinement  #
#--------------#
  0
#-------------------#
#  Copy boundaries  #
#-------------------#
  0
#-------------#
#  Smoothing  #
#-------------#
  1
    1  
    20  0.9
    0.0  0.0  0.0   1.0   1.0   1.0
