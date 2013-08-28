set terminal x11 size 1000, 500  #pngcairo solid enhanced font "arial,10" size 1000, 500 
cd 'pictures'
#set output '10000q0.png'
cd '..'
cd 'data'

lens  = "lens_shape/cylindrical_0.00625_0.00625_0.01.lns"
file1 = "1cylindrical_0.00625_0.00625_0.01_100lps_19.25A_50000eV_1000it_Cart.dat"

set xlabel "z (cm)"
set ylabel "r (mm)"

#set yrange [0:.28]
#set xrange [1:]

plot file1 using ($7):($6) title file1 with lines, \
     #file1 using ($3*100):($6) title file1 with lines, \

     
