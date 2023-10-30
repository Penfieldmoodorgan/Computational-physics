#----------------------------------------------------------------------------------

# Figura 1- Posicions pistons 1 i 5 en funció del temps

#----------------------------------------------------------------------------------
set term png
set output "P2-18P-fig1.png"
set key top left
set xlabel "temps (s)"
set ylabel "posició pistó(cm)"
plot "P2-18P-res1.dat" u 1:2 pt 7 title " P1",\
"P2-18P-res1.dat" u 1:6 pt 7 title " P5"

#----------------------------------------------------------------------------------

# Figura 2- Posicions pistons 3 i 5 en funció de la del pistó 2

#----------------------------------------------------------------------------------
set term png
set output "P2-18P-fig2.png"
set key top left
set xlabel "Posició pistó 2 (cm)"
set ylabel "Posició pistó(cm)"
plot "P2-18P-res1.dat" u 3:4 pt 7 title " P3",\
"P2-18P-res1.dat" u 3:6 pt 7 title " P5"

#-------------------------------------------------------------------------------------------------------

# Figura 3- Posició quart pistó calculada i interpolada, tant amb interpolació d'ordre zero com lineal

#-------------------------------------------------------------------------------------------------------
set term png
set output "P2-18P-fig3.png"
set key top right
set xlabel "t (s)"
set ylabel "Posició pistó 4 (cm)"
plot "P2-18P-res2.dat" u 1:2 w l title "int zero",\
"P2-18P-res2.dat" u 1:3 pt 3 ps 0 title "int lin",\
"P2-18P-res1.dat" u 1:5 pt 7 ps 1 title "X4 (t)"