set terminal png size 1366,768
set output 'BMEfunctionsStable.png'
set title "Porownanie y(t) z wartościami analitycznymi w warunkach numerycznej stabilności"

set xrange [-0.1 : 1.1]
set yrange [-0.5 : 1.5]
set xlabel "t"
set ylabel "y(t)"

set datafile separator whitespace

set grid

plot "BME.txt" using 1:2 every 10 title "y(t) approx" with points ls 1, \
     "BME.txt" using 1:3 title "analytical" with lines ls 2