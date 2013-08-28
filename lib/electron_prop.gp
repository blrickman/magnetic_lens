set terminal x11 size 1000, 500  #pngcairo solid enhanced font "arial,10" size 1000, 500 
cd 'pictures'
#set output '10000q0.png'
cd '..'
cd 'data'

lens  = "lens_shape/conical_0.01_0.0025_0.01.lns"
file1 = "e1-cylindrical_0.00625_0.00625_0.01_100lps_19.25A_50000eV_100000it_Cart.dat"
file2 = "e1-cylindrical_0.00625_0.00625_0.01_100lps_19.25A_50000eV_10000it_Cart.dat"

set xlabel "z (cm)"
set ylabel "r (mm)"

#set yrange [0:.28]
#set xrange [0:3]

plot lens  using ($1*100):($2*1000) notitle lc rgb "black" with lines,\
     file1 using ($3*100):(($4**2+$5**2)**.5*1000) title file1 with lines, \
     file1 using ($3*100):($6) title file2 with lines, \
