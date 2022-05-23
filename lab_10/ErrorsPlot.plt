set terminal png size 1366,768
set output 'Errors.png'
set title "wartości błędów bezwzględnych w zależności od kroku dt w metodzie "

set xrange [-10 : 0]
set yrange [-15 : 0]
set xlabel "dt"
set ylabel "|error|"

set datafile separator whitespace

set grid

plot "Errors.txt" using 1:2 title "BME" with linespoints ls 1, \
    "Errors.txt" using 1:3 title "PME" with linespoints ls 2, \
    "Errors.txt" using 1:4 title "PMT" with linespoints ls 1, x, 2*x