set terminal png size 1366,768
set output 'PMTfunctions.png'
set title "Porownanie y(t) z wartościami analitycznymi"

set xrange [-0.1 : 1.1]
set yrange [-0.5 : 1.5]
set xlabel "t"
set ylabel "y(t)"

set datafile separator whitespace

set grid

plot "PMT.txt" using 1:2 every 10 title "y(t) approx" with points ls 1, \
     "PMT.txt" using 1:3 title "analytical" with lines ls 2