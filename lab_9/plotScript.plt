#Wykres błedy
set xrange [-4.1 : -0.5]
set yrange [-15 : -3]
set terminal png size 1366,768
set output "Wykres_bledow.png"
set title 'Wykres bledow'
set ylabel 'log10(blad)'
set xlabel 'log10(krok)'

set grid
plot \
 "numErr.txt" using 1:2 with lines title "Bledy - dyskretyzacja Trojpunktowa", 4*x, 2*x