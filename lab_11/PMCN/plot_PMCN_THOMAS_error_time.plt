set xrange [-0.2 : 2.2]
 set yrange [ -7 : 0]
 set terminal png size 1366,768
 set output "error_time_PMCN_THOMAS.png"
 set title 'Zaleznosc logarytmu z błedu bezwględnego od czasu'
 set ylabel 'log10(error)'
 set xlabel 't'
 set grid
 plot \
   "PMCN_THOMAS_error_time.txt" using 1:2 every 15 with points title "log10(error)"