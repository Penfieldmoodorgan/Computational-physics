# gnuplot pre2-2
set term png
set output "P2-18P-fig2.png"
set key top left
set xlabel "Posició pistó 2 (cm)"
set ylabel "Posició pistó(cm)"
plot "P2-18P-res1.dat" u 3:4 pt 7 title " P3",\
"P2-18P-res1.dat" u 3:6 pt 7 title " P5"