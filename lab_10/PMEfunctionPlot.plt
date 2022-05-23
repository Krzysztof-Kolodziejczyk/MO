set terminal png size 1366,768
set output 'PMEfunctions.png'
set title "Porownanie y(t) z wartościami analitycznymi"

set xrange [-0.1 : 1.1]
set yrange [-0.5 : 1.5]
set xlabel "t"
set ylabel "y(t)"

set datafile separator whitespace

set grid

plot "PME.txt" using 1:2 title "y(t) approx" with linespoints ls 1, \
     "PME.txt" using 1:3 title "analytical" with lines ls 2