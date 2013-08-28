set terminal x11 size 1000, 500  #pngcairo solid enhanced font "arial,10" size 1000, 500 
cd 'pictures'
#set output '10000q0.png'
cd '..'
cd 'data'

fn = "cylindrical_0.00625_0.00625_0.01_100lps_19.25A_50000eV_1000it_Cart"

lens  = "lens_shape/" . "cylindrical_0.00625_0.00625_0.01.lns"
file1 = "e1-" . fn . ".dat"
file2 = "e2-" . fn . ".dat"
file3 = "e3-" . fn . ".dat"

beam1(x) = w01 * sqrt(1 + (x - xc1)**2 / xR1**2 )
beam2(x) = w02 * sqrt(1 + (x - xc2)**2 / xR2**2 )
beam3(x) = w03 * sqrt(1 + (x - xc3)**2 / xR3**2 )

min = 430
max = 440
skip = 30

set fit quiet
FIT_LIMIT = 1e-6
set fit logfile "fit_log/" . fn

plot file1 using ($3*100):(($1**2+$2**2)**.5*1000) every ::min::max
pause -1 "Hit return to continue"

fit beam1(x) file1 using ($3*100):(($1**2+$2**2)**.5*1000) every ::min::max via w01,xc1,xR1
fit beam2(x) file2 using ($3*100):(($1**2+$2**2)**.5*1000) every ::min::max via w02,xc2,xR2
fit beam3(x) file3 using ($3*100):(($1**2+$2**2)**.5*1000) every ::min::max via w03,xc3,xR3

print "Fit Results:"
print "Beam Waist for file1 is " , w01 , " and focus occurs at " , xc1
print "Beam Waist for file2 is " , w02 , " and focus occurs at " , xc2
print "Beam Waist for file3 is " , w03 , " and focus occurs at " , xc3

set xlabel "z (cm)"
set ylabel "r (mm)"

set yrange [0:.28]
set xrange [0.5:3]

plot lens  using ($1*100):($2*1000) notitle lc rgb "black" with lines,\
     file1 using ($3*100):(($1**2+$2**2)**.5*1000) every skip title 'file1' with points, \
     file2 using ($3*100):(($1**2+$2**2)**.5*1000) every skip title 'file2' with points, \
     file3 using ($3*100):(($1**2+$2**2)**.5*1000) every skip title 'file3' with points, \
     beam1(x) title "Best Fit for file1" with lines, \
     beam2(x) title "Best Fit for file2" with lines, \
     beam3(x) title "Best Fit for file3" with lines
