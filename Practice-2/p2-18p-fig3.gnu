#graf 3 prac 2
set term png
set output "P2-18P-fig3.png"
set key top right
set xlabel "t (s)"
set ylabel "Posició pistó 4 (cm)"
plot "P2-18P-res2.dat" u 1:2 w l title "int zero",\
"P2-18P-res2.dat" u 1:3 pt 3 ps 0 title "int lin",\
"P2-18P-res1.dat" u 1:5 pt 7 ps 1 title "X4 (t)"