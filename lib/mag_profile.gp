set terminal x11 size 1000, 500  #pngcairo solid enhanced font "arial,10" size 1000, 500 

file = 'conical_0.01_0.0025_99_0.01_10_200'
cd 'solenoid-data/' . file

set xlabel "z (cm)"
set ylabel "B (T)"

#unset key

#set yrange [0:.255]
#set xrange [6.07:6.09]

plot for [fn in system("ls *.dat")] fn using ($2*100):($3) title 'magnetic field along ' . fn with lines
