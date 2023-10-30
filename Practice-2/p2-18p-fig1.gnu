# gnuplot pr2-1
set term png
set output "P2-18P-fig1.png"
set key top left
set xlabel "temps (s)"
set ylabel "posició pistó(cm)"
plot "P2-18P-res1.dat" u 1:2 pt 7 title " P1",\
"P2-18P-res1.dat" u 1:6 pt 7 title " P5",\