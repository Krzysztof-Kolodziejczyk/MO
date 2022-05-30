set terminal png size 1366,768
set output 'functions.png'
set title "Porownanie U(x) wartościami analitycznymi z przyblizonymi"

set xrange [-0.1 : 1.1]
set yrange [-0.5 : 1.1]
set xlabel "x"
set ylabel "U(x)"

set datafile separator whitespace

set grid

plot "numFun.txt" using 1:2 title "U(x) analytical" with lines, \
     "numFun.txt" using 1:3 every 3 title "Ux approx numerow" with points ls 4,\
     "conFun.txt" using 1:2 every 3 title "Ux approx conventional" with points ls 3,\
