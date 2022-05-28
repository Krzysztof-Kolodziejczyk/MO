set xrange [-2 : -0.2]
set yrange [ -6 : -2]
set terminal png size 1366,768
set output "error_PMCN_THOMAS.png"
set title 'Zaleznosc logarytmu z błedu bezwględnego od logarytmu z korku h'
set ylabel 'log10(error)'
set xlabel 'log10(h)'
set grid
plot \
 "PMCN_THOMAS_error.txt" using 1:2 with points pt 1 lc 10 title "KMB",\
  "PMCN__THOMAS_error.txt" using 1:2 with lines title "KMB", 2 *x