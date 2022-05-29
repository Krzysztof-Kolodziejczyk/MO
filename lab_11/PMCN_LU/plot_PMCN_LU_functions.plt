set xrange [0 : 12]
set yrange [ -0.2 : 1.2]
set terminal png size 1366,768
set output "functions_PMCN_LU.png"
set title 'rozwiązania numeryczne i analityczne w chwilach czasowych t = 0.5, t = 1.0, t = 2.0 dla PMCN z zastosowaniem dekompozycji LU'
set ylabel 'U(x,t)'
set xlabel 'x'
set grid
plot \
  "PMCN_LU_functions.txt" using 1:2 with lines title "analytical t = 0.5",\
  "PMCN_LU_functions.txt" using 1:2 every 15 with points pt 1 lc 9 title "KMB t = 0.5",\
  "PMCN_LU_functions.txt" using 1:4 with lines title "analytical t = 1.0",\
  "PMCN_LU_functions.txt" using 1:5 every 15 with points pt 1 lc 6 title "KMB t = 1.0",\
  "PMCN_LU_functions.txt" using 1:6 with lines title "analytical t = 2.0",\
  "PMCN_LU_functions.txt" using 1:7 every 15 with points pt 1 lc 8 title "KMB t = 2.0"