set terminal x11 size 1000, 500  #pngcairo solid enhanced font "arial,10" size 1000, 500 
cd 'pictures'
#set output '10000q0.png'
cd '../data'

dir = "cyl-f_cyl-b_0.00625_0.00625_0.01_100lps_18.9995A_s2-0.01904848498_50000eV_100000it_Cart"

lens1  = "lens_shape/" . "conical-front_0.01_0.0025_0.01.lns"
lens2  = "lens_shape/" . "conical-back_0.0025_0.01_0.01.lns"
file1 = dir . "/e1.dat"
file2 = dir . "/e2.dat"
file3 = dir . "/e3.dat"

beam1(x) = w01 * sqrt(1 + (x - xc1)**2 / xR1**2 )
beam2(x) = w02 * sqrt(1 + (x - xc2)**2 / xR2**2 )
beam3(x) = w03 * sqrt(1 + (x - xc3)**2 / xR3**2 )

linear1(x) = m1 * x + b1
linear2(x) = m2 * x + b2
linear3(x) = m3 * x + b3

scalx = 100
scaly = 1000
unitx = 'cm'
unity = 'mm'
min = 42500
max = 43500
lmin = 78000
skip = 30
lskip = 3000

beamslope(x,w0,xc,xR) = w0 * (x - xc) / xR**2 / sqrt(1 + (x - xc)**2 / xR**2) * scalx / scaly

set fit quiet
FIT_LIMIT = 1e-6
logfn = "fit_log/" . dir . ".log"
resultslog = "fit_log/" . dir . "-RESULTS.log"
set fit logfile logfn

#if (strlen(system("file " . logfn . " | grep ERROR")) == 0 ) {
#  load logfn
#  print m1
#}

plot file1 using ($3*scalx):(($1**2+$2**2)**.5*scaly) every skip::min::max
#pause -1 "Hit return to continue"

plot file1 using ($3*scalx):(($1**2+$2**2)**.5*scaly) every lskip::lmin
#pause -1 "Hit return to continue"

fit beam1(x) file1 using ($3*scalx):(($1**2+$2**2)**.5*scaly) every ::min::max via w01,xc1,xR1
print "Fit Results:"
print "The beam waist for file1 is " , sprintf('%.4f',w01/scaly/10**-6) , "um and focus occurs at " , sprintf('%.5f',xc1/scalx/10**-2) , "cm with a Ralaigh range of " , sprintf('%.5f',xR1/scalx/10**-6) , "um"

w02 = w01; xc2 = xc1; xR2 = xR1;
fit beam2(x) file2 using ($3*scalx):(($1**2+$2**2)**.5*scaly) every ::min::max via w02,xc2,xR2
print "The beam waist for file2 is " , sprintf('%.4f',w02/scaly/10**-6) , "um and focus occurs at " , sprintf('%.5f',xc2/scalx/10**-2) , "cm with a Ralaigh range of " , sprintf('%.5f',xR2/scalx/10**-6) , "um"

w03 = w01; xc3 = xc1; xR3 = xR1;
fit beam3(x) file3 using ($3*scalx):(($1**2+$2**2)**.5*scaly) every ::min::max via w03,xc3,xR3
print "The beam waist for file3 is " , sprintf('%.4f',w03/scaly/10**-6) , "um and focus occurs at " , sprintf('%.5f',xc3/scalx/10**-2) , "cm with a Ralaigh range of " , sprintf('%.5f',xR3/scalx/10**-6) , "um"

print ""
print "Aberrations: " , sprintf('%.5f',(xc1-xc2)/scalx/10**-6) , "um and " , sprintf('%.5f',(xc2-xc3)/scalx/10**-6) ,"um"
print "Beam Slopes: " , sprintf('%.4f',1000*beamslope(.012*scalx,w01,xc1,xR1)), "mrad, " , sprintf('%.4f',1000*beamslope(.012*scalx,w02,xc2,xR2)), "mrad, " , sprintf('%.4f',1000*beamslope(.012*scalx,w03,xc3,xR3)) , "mrad"

fit linear1(x) file1 using ($3*scalx):(($1**2+$2**2)**.5*scaly) every 100*skip::lmin via m1,b1
fit linear2(x) file2 using ($3*scalx):(($1**2+$2**2)**.5*scaly) every 100*skip::lmin via m2,b2
fit linear3(x) file3 using ($3*scalx):(($1**2+$2**2)**.5*scaly) every 100*skip::lmin via m3,b3

print "Second Foci: " , sprintf('%.3f',-b1/m1/scalx), "m, " , sprintf('%.3f',-b2/m2/scalx), "m, "  , sprintf('%.3f',-b3/m3/scalx), "m"

set xlabel "z (". unitx .")"
set ylabel "r (". unity .")"

set yrange [0:.28]
set xrange [.005*scalx:.025*scalx]

plot lens1  using ($1*scalx):($2*scaly) notitle lc rgb "black" with lines,\
     lens2  using ($1*scalx):($2*scaly) notitle lc rgb "black" with lines,\
     file1 using ($3*scalx):(($1**2+$2**2)**.5*scaly) every skip title 'file1' with points, \
     file2 using ($3*scalx):(($1**2+$2**2)**.5*scaly) every skip title 'file2' with points, \
     file3 using ($3*scalx):(($1**2+$2**2)**.5*scaly) every skip title 'file3' with points, \
     beam1(x) title "Best Fit for file1" with lines, \
     beam2(x) title "Best Fit for file2" with lines, \
     beam3(x) title "Best Fit for file3" with lines, \
     linear1(x) title "linear fit" with lines, \
     linear2(x) title "linear fit" with lines, \
     linear3(x) title "linear fit" with lines

set print resultslog
print strftime("%H:%M:%.3S %d-%b-%Y",time(0.0))
print ""
print "Focus fit (w0,xc,xR) for:"
print "Beam 1: " , w01, xc1, xR1
print "Beam 2: " , w02, xc2, xR2
print "Beam 3: " , w03, xc3, xR3
print ""
print "Linear fit (m,b) for:"
print "Beam 1: " , m1, b1
print "Beam 2: " , m2, b2
print "Beam 3: " , m3, b3
print ""
print "Aberrations: " , sprintf('%.5f',(xc1-xc2)/scalx/10**-6) , "um and " , sprintf('%.5f',(xc2-xc3)/scalx/10**-6) ,"um"
print "Beam Slopes: " , sprintf('%.4f',1000*beamslope(.012*scalx,w01,xc1,xR1)), "mrad, " , sprintf('%.4f',1000*beamslope(.012*scalx,w02,xc2,xR2)), "mrad, " , sprintf('%.4f',1000*beamslope(.012*scalx,w03,xc3,xR3)) , "mrad"
print "Second Foci: " , sprintf('%.3f',-b1/m1/scalx), "m, " , sprintf('%.3f',-b2/m2/scalx), "m, "  , sprintf('%.3f',-b3/m3/scalx), "m"
