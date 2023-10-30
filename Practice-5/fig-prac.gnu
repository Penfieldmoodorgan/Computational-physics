#----------------------------------------------------------------------------------

# Figura 1- Histograma normalitzat amb 120 caixes comparat amb la distribució exacta

#----------------------------------------------------------------------------------
set term png
set output "P5-18P-fig1.png"
set key top right
set xlabel "x"
set ylabel "P(x)"
set grid
plot "P5-18P-res1.dat" index 0 u 1:2:3 w yerrorbars lc rgb "#C51D34" notitle,\
"P5-18P-res1.dat" index 0 u 1:2 w p pt 7 lc rgb "#C51D34" notitle,\
"P5-18P-res1.dat" index 0 u 1:2 with histeps lc rgb "#C51D34" notitle,\
"P5-18P-res1.dat" index 1 u 1:2 w l lc rgb "#641C34" title "Exacte"

#----------------------------------------------------------------------------------

# Figura 2-Trajectòria durant tota l'evolució de les 5 molècules

#----------------------------------------------------------------------------------
set term png
set output "P5-18P-fig2.png"
set key top left
set xlabel "x(m)"
set ylabel "y(m)"
set grid
plot "P5-18P-res2.dat" index 0 u 2:3 w l title "Particula 1",\
"P5-18P-res2.dat" index 0 u 4:5 w l title "Particula 2",\
"P5-18P-res2.dat" index 0 u 6:7 w l lc rgb "red" title "Particula 3",\
"P5-18P-res2.dat" index 0 u 8:9 w l title "Particula 4",\
"P5-18P-res2.dat" index 0 u 10:11 w l lc rgb "green" title "Particula 5"

#----------------------------------------------------------------------------------

# Figura 3-t,Var(y(t)) comparat amb 2Dt

#----------------------------------------------------------------------------------
set term png
set output "P5-18P-fig3.png"
set key top left
set xlabel "t"
set ylabel "Var(y)"
set grid
plot "P5-18P-res2.dat" index 1 u 1:2 w l title "Var(y(t))"



