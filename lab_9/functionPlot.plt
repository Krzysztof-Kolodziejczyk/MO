set terminal png size 1366,768
set output 'functions.png'
set title "Porownanie U(x) wartościami analitycznymi z przyblizonymi"

set xrange [-0.1 : 1.1]
set yrange [-0.5 : 1.1]
set xlabel "x"
set ylabel "U(x)"

set datafile separator whitespace

set grid

plot "numFun.txt" using 1:2 title "U(x) true" with linespoints ls 1, \
     "numFun.txt" using 1:3 title "Ux approx" with linespoints ls 2