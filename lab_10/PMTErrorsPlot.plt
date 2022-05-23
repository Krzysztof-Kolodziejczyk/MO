set terminal png size 1366,768
set output 'PMTErrors.png'
set title "wartości błędów bezwzględnych w zależności od kroku dt w metodzie PMT"

set xrange [-10 : 0]
set yrange [-15 : 0]
set xlabel "dt"
set ylabel "|error|"

set datafile separator whitespace

set grid

plot "PMTErrors.txt" using 1:2 title "error" with linespoints ls 1, x