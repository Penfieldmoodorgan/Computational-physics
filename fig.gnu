#----------------------------------------------------------------------------------

# Figura 1- Sumatori en funció de N i comportament asimptòtic

#----------------------------------------------------------------------------------
set term png
set output "P1-18P-fig1.png"
set key bottom right
set xlabel "N"
set ylabel "SUMA(8,N)"
set logscale y
plot "P1-18P-res1.dat" u 1:2 pt 7 title "SUMA(8,N) (N)",\
"P1-18P-res1.dat" u 1:2 w l title "comportament asimptòtic"

#----------------------------------------------------------------------------------

# Figura 2- Sumatori entre comportament asimptòtic funció de N 

#----------------------------------------------------------------------------------
set term png
set output "P1-18P-fig2.png"
set key top right
set xlabel "N"
set ylabel "SUMA(8,N)/SUMA asimpt"
plot "P1-18P-res2.dat" u 1:2 pt 7 title "SUMA(8,N)/SUMA asimpt (N)"
