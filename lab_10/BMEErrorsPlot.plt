set terminal png size 1366,768
set output 'BMEErrors.png'
set title "wartościami błędów bezwzględnych w zależności od kroku dt"

set xrange [-0.1 : 1.1]
set yrange [-0.5 : 1.5]
set xlabel "dt"
set ylabel "|error|"

set datafile separator whitespace

set grid

plot "BMEErrors.txt" using 1:2 title "y(t) approx" with linespoints ls 1,