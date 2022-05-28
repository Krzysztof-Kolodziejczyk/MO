set xrange [-2 : -0.2]
set yrange [ -6 : -2]
set terminal png size 1366,768
set output "error_PMCN_LU.png"
set title 'Zaleznosc logarytmu z błedu bezwględnego od logarytmu z korku h'
set ylabel 'log10(error)'
set xlabel 'log10(h)'
set grid
plot \
 "PMCN_LU_error.txt" using 1:2 with points pt 1 lc 10 title "PMCN LU",\
  "PMCN_LU_error.txt" using 1:2 with lines title "PMCN_THOMAS LU", 2 *x