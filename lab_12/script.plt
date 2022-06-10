set xrange [-1.2 : 1.2]
set yrange [ -0.2 : 1.2]
set terminal png size 1366,768
set output "function.png"
set title 'funkcja'
set ylabel 'f(x)'
set xlabel 'x'
set grid
plot \
  "function.txt" using 1:2 with lines title "function",\
  "normal.txt" using 1:2 with points pt 1 lc 10 title "normal",\
  "czebyszew.txt" using 1:2 with points pt 2 lc 8 title "czebyszew",\