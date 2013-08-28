set terminal x11 size 1000, 1000  #pngcairo solid enhanced font "arial,10" size 1000, 500 
cd 'pictures'
#set output '10000q0.png'
cd '..'
cd 'data'

file2 = "conical_0.01_0.0025_0.01_100lps_20A_50000eV_10000it_Cart.dat"
set view equal xy
set xlabel "x (mm)"
set ylabel "y (mm)"
set zlabel "z (cm)"

set xtics nomirror
set ztics -5,2.5,10
set ytics -2,1,2
set xtics -2,1,2

#set zrange [-5:10]
set xrange [-2:2]
set yrange [-2:2]

splot file2 using ($1*1000):($2*1000):($3*100) title file2 with lines, \
