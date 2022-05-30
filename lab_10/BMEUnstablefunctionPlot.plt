set terminal png size 1366,768
set output 'BMEfunctionsUnstable.png'
set title "Porownanie y(t) z wartościami analitycznymi w warunkach numerycznej niestabilności"

set xrange [-0.1 : 1.1]
set yrange [-30 : 30]
set xlabel "t"
set ylabel "y(t)"

set datafile separator whitespace

set grid

plot "BME_unstable.txt" using 1:2 title "y(t) approx" with points ls 1, \
     "BME_unstable.txt" using 1:3 title "analytical" with lines ls 2