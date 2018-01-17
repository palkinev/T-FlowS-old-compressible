#!/bin/bash

# this script generated test cases for Problem1, divide domain in div_subdomains parts,
# launches TFlowS, extract final precision and plot it in xmgrace

# Generate, (Divide) must be compiled

#-----------------------------

div_subdomains=1
div_method=INE
max_power=6
xmgrace=0
octave=1

#-----------------------------
Tdir="./"
cp /dev/null L2_Error.dat # clean log

for (( i = 6; i <= $max_power; i++ ))
do

	Tdir="$(echo "2^$i" | bc -l)"/

	if [ ! -d "$Tdir" ]; then # if dir doesn't exist.
    	mkdir "$Tdir"
    fi

	cd "$Tdir"

#-----------------------------

#create gen_in file
cat << EOF > gen_in
gen
0
0
x
yz
0
0
0
1 1 4
C
BC
skip

EOF

#-----------------------------

#create gen.d file
cat << EOF > gen.d
#
#     8-----------6
#    /|          /|
#   7-----------5 |
#   | |         | |
#   | |         | |
#   | |         | |
#   | 4 - - - - | 2
#   |/          |/
#   3-----------1
#
#   x: 1 -> 2
#   y: 1 -> 3
#   z: 1 -> 5
#
#-------------------------------------------#
#  Nodes (cells), boundary cells and sides  #
#-------------------------------------------#
  30000 10000 90000

#----------#
#  Points  #
#----------#
8
  1   0.0            -3.0            -3.0
  2   2.0            -3.0            -3.0
  3   0.0             3.0            -3.0
  4   2.0             3.0            -3.0
  5   0.0            -3.0             3.0
  6   2.0            -3.0             3.0
  7   0.0             3.0             3.0
  8   2.0             3.0             3.0        
#----------#
#  Blocks  #
#----------#
1
  1 $(echo "2^$i+1" | bc -l)  4  4 
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
  2
    1     Imin 
        1   1
    2     Imax
        1   2


#-----------------------#
#  Periodic boundaries  #
#-----------------------#
  2
      3    1 2 4 3
           5 6 8 7
      4    1 2 6 5
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

EOF

../../../Generate/gen < gen_in

#-----------------------------

#create div_in file
cat << EOF > div_in
gen
$div_subdomains
$div_method
skip
EOF

if [ $div_subdomains != 1 ]; then
  ../../../Divide/divide < div_in
fi

#-----------------------------

#create gen.b file
cat << EOF > gen.b
  1            #VISc
  1   FLUID    0.03  1.00   1.0    1.0
  2
   1       INFLOW       0.0   0.0     0.0        0.0   0.0     0.0   0.0   0.0   0.0
   2       OUTFLOW      0.0   0.0     0.0        0.0   0.0     0.0   0.0   0.0   0.0   0.0     0.0   0.0
  1
  1   0.0   0.0     0.0   0.0  0.0

EOF

#-----------------------------

#create T-Flows.cmn
cat << EOF > T-FlowS.cmn
gen
MOIN_1_T HOT
skip
800
189000
1
0.1 0.1 0.1
0.1 0.1 0.1
DNS
NO
SIMPLE
1.0 # U    % URF
1.0 # P    % URF
1.0 # Zmix % URF
par
fi
fi
fi
CDS
CDS
ic
1.e-8 # U % STol
1.e-8 # P % STol
1.e-8 # Zmix % STol
1.e-7 # SIMTol
0.00125
0.0 0.0 0.0
0.0 0.0 0.0
skip
skip
skip
skip
skip
EOF

#-----------------------------

# copy stdout to
exec > >(tee -i output.log)
# also, append stderr to it
exec 2>&1

if [ $div_subdomains == 1 ]; then
  ../../../Process/TFLowS
else
  mpirun.mpich -np 4 ../../../Process/TFLowS
fi

#-----------------------------

# filter output.log
awk '/Error/{nr[NR+2]}; NR in nr' output.log > L2.log
sed -i -e "s%#\ %%g" L2.log

cd ../

#rm -r "$Tdir"
    
done

#-----------------------------

if [ $xmgrace == 1 ]; then
# create tmp.dat
cat << EOF > tmp.dat
1 0
1 1
EOF

#create grace.bat
cat << EOF > grace.bat
AUTOSCALE ONREAD NONE

read block "L2_Error.dat"

#U
block xy "1:3"
s0 LEGEND "U"
s0 LINE TYPE 0
s0 SYMBOL 1
s0 SYMBOL SIZE 0.98
s0 SYMBOL COLOR 2
s0 SYMBOL FILL 1
s0 SYMBOL FILL COLOR 2
#V
block xy "1:4"
s1 LEGEND "V"
s1 LINE TYPE 0
s1 SYMBOL 2
s1 SYMBOL SIZE 0.98
s1 SYMBOL COLOR 3
s1 SYMBOL FILL 1
s1 SYMBOL FILL COLOR 3
#P
block xy "1:5"
s2 LEGEND "P"
s2 LINE TYPE 0
s2 SYMBOL 2
s2 SYMBOL SIZE 0.98
s2 SYMBOL COLOR 4
s2 SYMBOL FILL 1
s2 SYMBOL FILL COLOR 4
#Z
block xy "1:6"
s3 LEGEND "Z"
s3 LINE TYPE 0
s3 SYMBOL 2
s3 SYMBOL SIZE 0.98
s3 SYMBOL COLOR 5
s3 SYMBOL FILL 1
s3 SYMBOL FILL COLOR 5
#Rho
block xy "1:7"
s4 LEGEND "Rho"
s4 LINE TYPE 0
s4 SYMBOL 2
s4 SYMBOL SIZE 0.98
s4 SYMBOL COLOR 6
s4 SYMBOL FILL 1
s4 SYMBOL FILL COLOR 6

#1nd order line
read block "tmp.dat"
block xy "1:2"
kill block
s5.x[0] = s2.x[imin(s2.x)]
s5.x[1] = s2.x[imax(s2.x)]

# based on pressure
#s5.y[0] = s2.y[0]
# based on  rho
s5.y[0] = s4.y[0]
s5.y[1] = s5.y[0]/(s2.x[imax(s2.x)]/s2.x[imax(s2.x)-1])

#2nd order line
copy s5 to s6
s6.y[1] = s6.y[0]/(s2.x[imax(s2.x)]/s2.x[imax(s2.x)-1])^2


# appearence
world xmin 32*0.99
world xmax 64*1.01
world ymin s0.y[imax(s2.x)]*0.99
world ymax s5.y[0]*1.01
yaxes scale logarithmic
yaxis  ticklabel format exponential
yaxis  ticklabel prec 0


xaxis  tick major 32
yaxis  tick major 10
xaxis  ticklabel char size 2.000000
yaxis  ticklabel char size 2.000000
legend 1.689, 0.99

s5 LEGEND "1st order"
s5 line type 1
s5 linewidth 3
s5 linestyle 3
s5 line color 1

s6 LEGEND "2nd order"
s6 line type 1
s6 linewidth 3
s6 linestyle 3
s6 line color 1


# print
PRINT TO "L2_Error.png"
HARDCOPY DEVICE "PNG"
DEVICE "PNG" DPI 300
PAGE SIZE 1980, 1050
#XMIN YMIN XMAX YMAX
VIEW 0.14, 0.09, 1.88, 0.99
DEVICE "PNG" FONT ANTIALIASING on
PRINT

SAVEALL "L2_Error.agr"

EOF

gracebat -batch grace.bat -nosafe
fi

#-----------------------------

if [ $octave == 1 ]; then
cat << EOF > tmp.m
close all hidden
graphics_Toolkit('gnuplot'); % set gnuplot as plot program

% see also: available_graphics_toolkits

% A(i,j)
% i = test id
% j = 1 -> mesh nx
% j = 2 -> U tolerance
% j = 3 -> V tolerance
% j = 4 -> P tolerance
% j = 5 -> Z tolerance
% j = 6 -> Rho tolerance
fid = fopen ("L2_Error.dat", "r");
A = transpose(fscanf(fid,"%e",[7 Inf]));
fclose(fid);
A(:,2) = [];

% first order accurate
first = A(size(A(:,1))+find(A(:,2:end)==max(max(A(:,2:end)))))*1.01;
first(2) = first(1)/2;
% second order accurate
second = first;
second(2) = second(1)/4;


% plot
h = semilogy(...
  A(:,1),A(:,2),'ro','MarkerFaceColor','r',...
  A(:,1),A(:,3),'gd','MarkerFaceColor','g',...
  A(:,1),A(:,4),'bp','MarkerFaceColor','b',...
  A(:,1),A(:,5),'mv','MarkerFaceColor','m',...
  A(:,1),A(:,6),'co','MarkerFaceColor','c',...
  A([1 end],1),first,'k--','LineWidth',4,...
  A([1 end],1),second,'k--','LineWidth',4);  
set(h, 'MarkerSize', 20);

% appearence
H = 21; W = 29.7; % A4 paper size
set(gcf(),'PaperUnits','centimeters')
set(gcf(),'PaperSize',[H,W])
set(gcf(),'PaperPosition',[0,0,W,H])

%change axis properties
h = get (gcf, "currentaxes");
set(h,"fontweight","bold")
set(h, 'FontSize', 14);
title('Problem 1 test case');
xlabel('mesh Nx');
ylabel('precision');
set(h,'XTick',A(:,1));
set(h,'XtickLabel',A(:,1));
% set(gca,'YtickLabel',sprintf('%.0e|',A(:,2)))

%change legend properties
h = legend('U','V','P','Z','Rho','Location','southwest');
set(h,'FontSize',10);
set(h,'PlotBoxAspectRatioMode','manual');
set(h,'PlotBoxAspectRatio',[1 0.00000000000001 1]);

%export to file
%print("L2_Error.eps", "-deps","-r300","-S1920,1080","-color");
print("L2_Error.eps", "-deps","-color");

EOF

octave --no-gui tmp.m

# remove all white
#convert -density 300 -transparent white L2_Error.eps L2_Error.png

# remove all transparent
convert -density 300 -flatten L2_Error.eps L2_Error.png  
fi