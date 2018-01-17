#!/bin/bash

# this script generated test cases for Problem1, divide domain in div_subdomains parts,
# launches TFlowS, extract final precision and plot it in xmgrace

# Generate, (Divide) must be compiled

#-----------------------------

div_subdomains=1
div_method=INE
min_power=6 # Moin value - 6
max_power=6 # Moin value - 12, but max_power >=11 diverges with any Time_Steps
declare -a N_table
Time_Steps=(800) # Moin value (2 5 10 25 50 100 200 400 800)
clean_log=false
octave=true
octave=false


#-----------------------------
for Nt in ${Time_Steps[@]}
do


#-----------------------------
Tdir="./"

if [ ! -f L2_Error_$Nt.dat ]; then
    touch L2_Error_$Nt.dat
elif [ "$clean_log"=true ]; then
    cp /dev/null L2_Error_$Nt.dat # clean log
fi
sed -i -e '/\#/d' L2_Error_$Nt.dat

for (( i = $min_power; i <= $max_power; i++ ))
do

Tdir="$(echo "2^$i" | bc -l)"/

if [ ! -d "$Tdir" ]; then # if dir doesn't exist.
    mkdir "$Tdir"
fi

cd "$Tdir"

#-----------------------------

#if [ $Nt == 800 ] || [ $i == 11 ]; then
  
echo "Nx = ", "$(echo "2^$i" | bc -l)", "Nt = ", $Nt

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
  3000000 1000000 9000000

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
  1
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
MOIN_1 HOT
skip
$(echo "$Nt" | bc -l)
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
1.e-9 # U % STol
1.e-9 # P % STol
1.e-9 # Zmix % STol
1.e-8 # SIMTol
$(echo "1/$Nt" | bc -l)
0.0 0.0 0.0
0.0 0.0 0.0
skip
skip
skip
skip
skip
EOF
#-----------------------------

# copy stdout to output.log
#exec &> >(tee -ia output.log)
# also, append stderr to it
#exec 2>&1

if [ $div_subdomains == 1 ]; then
  ../../../Process/TFLowS > output.log
else
  mpirun.mpich -np $div_subdomains ../../../Process/TFLowS > output.log
fi

#-----------------------------

# filter output.log
awk '/Error/{nr[NR+2]}; NR in nr' output.log > L2.log
sed -e "s%#\ %%g" L2.log > tmp.log
n=$(echo "2^$i" | bc -l)
awk -v prefix="$n   " '{print prefix $0}' tmp.log > L2.log

tail -n1 L2.log >> ../L2_Error_$Nt.dat

#fi

cd ../


done

#-----------------------------

done

#cat *.dat > L2_Error.dat
#cat << EOF > tmp.m
#close all hidden
#graphics_toolkit('gnuplot'); % set gnuplot as plot program
#
#% see also: available_graphics_toolkits
#
#% A(i,j)
#% i = test id
#% j = 1 -> mesh nx
#% j = 2 -> U tolerance
#% j = 3 -> V tolerance
#% j = 4 -> P tolerance
#% j = 5 -> Z tolerance
#% j = 6 -> Rho tolerance
#filename = "L2_Error.dat";
#
#if exist(filename,"file")
#  fid = fopen (filename, "r");
#  A = transpose(fscanf(fid,"%e",[7 Inf]));
#  fclose(fid);
#  A(:,2) = [];
#
#  % first order accurate
#  first = A(size(A(:,1))+find(A(:,2:end)==max(max(A(:,2:end)))))*1.01;
#  first(2) = first(1)/2;
#  % second order accurate
#  second = first;
#  second(2) = second(1)/4;
#
#
#  % plot
#  h = semilogy(...
#    A(:,1),A(:,2),'ro','MarkerFaceColor','r',...
#    A(:,1),A(:,3),'gd','MarkerFaceColor','g',...
#    A(:,1),A(:,4),'bp','MarkerFaceColor','b',...
#    A(:,1),A(:,5),'mv','MarkerFaceColor','m',...
#    A(:,1),A(:,6),'co','MarkerFaceColor','c',...
#    A([1 end],1),first,'k--','LineWidth',4,...
#    A([1 end],1),second,'k--','LineWidth',4);  
#  set(h, 'MarkerSize', 20);
#
#  % appearence
#  H = 21; W = 29.7; % A4 paper size
#  set(gcf(),'PaperUnits','centimeters')
#  set(gcf(),'PaperSize',[H,W])
#  set(gcf(),'PaperPosition',[0,0,W,H])
#
#  %change axis properties
#  h = get (gcf, "currentaxes");
#  set(h,"fontweight","bold")
#  set(h, 'FontSize', 14);
#  title('Problem 1 test case');
#  xlabel('mesh Nx');
#  ylabel('precision');
#  set(h,'XTick',A(:,1));
#  set(h,'XtickLabel',A(:,1));
#  % set(gca,'YtickLabel',sprintf('%.0e|',A(:,2)))
#
#  %change legend properties
#  h = legend('U','V','P','Z','Rho','Location','southwest');
#  set(h,'FontSize',10);
#  set(h,'PlotBoxAspectRatioMode','manual');
#  set(h,'PlotBoxAspectRatio',[1 0.00000000000001 1]);
#
#  %export to file
#  %print("L2_Error.eps", "-deps","-r300","-S1920,1080","-color");
#  print("L2_Error.eps", "-deps","-color");
#else
#  disp(strcat("file ", filename," does not exist"))
#endif
#
#
#EOF
#
#if [ "$octave"=true ] && [ -f L2_Error.eps ]; then
#    octave --no-gui tmp.m
#
#    # remove all white
#    #convert -density 300 -transparent white L2_Error.eps L2_Error.png
#
#    # remove all transparent
#    convert -density 300 -flatten L2_Error.eps L2_Error.png
#fi
