set terminal pngcairo solid enhanced font "arial,10" size 1000, 500 

set output '10000d0.png'

file1 = "/data/concical_1_025_1_100_10A-100000it.csv"

set xlabel "z (cm)"
set ylabel "r (mm)"

set xrange [6.07:6.09]

plot file1 using ($3*100):($1*1000) title 'electron path 100000 it' with lines#, \

show output

