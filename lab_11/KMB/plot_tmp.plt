set xrange [0 : 12]
set yrange [ -16 : 9]
set terminal png size 1366,768
set output "tmp.png"
set title 'Zaleznosc logarytmu z błedu bezwględnego od x'
set ylabel 'log10(error)'
set xlabel 'x'
set grid
plot \
 "tmp.txt" using 1:2 with points pt 1 lc 10 title "KMB"