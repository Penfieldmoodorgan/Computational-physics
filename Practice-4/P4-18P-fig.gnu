#----------------------------------------------------------------------------------

# Figura 1- Isoterma que conté només punts (V,P) estables (P'(V)<0)

#----------------------------------------------------------------------------------
set term png
set output "P4-18P-fig1.png"
set key top right
set xlabel "v(unitats reduides)"
set ylabel "P,dP(unitats reduides)"
plot "P4-18P-res.dat" index 0 u 1:($2<0?$3:1/0) w l title "P(v)"

#----------------------------------------------------------------------------------

# Figura 2- Corba a T=0.92 amb V€[1/3+0.1,2]

#----------------------------------------------------------------------------------
set term png
set output "P4-18P-fig2.png"
set key top right
set xlabel "v(unitats reduides)"
set ylabel "P(unitats reduides)"
plot "P4-18P-res.dat" index 1 u 1:2 w l lc rgb "#CF3476" title "P(v)",\
0 notitle lc black