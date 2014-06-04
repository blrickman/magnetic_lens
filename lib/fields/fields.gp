set terminal pngcairo solid font "times,18" size 1000, 750 

set datafile separator ","
set key title "Height in r"

set label 1 system("head -n5 " . datafile) at graph 0,1

set xlabel "z"
set ylabel "Field"
set yrange [-rangey:rangey]

plot for [i=2:height] datafile using 1:i w l title columnhead(i)
