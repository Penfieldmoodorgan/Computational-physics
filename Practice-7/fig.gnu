#----------------------------------------------------------------------------------

# Figura 1- Oscil·lacions petites

#----------------------------------------------------------------------------------
set term png
set output "P7-18P-b-fig1.png"
set key bottom box outside horizontal
set title "OSCIL·LACIONS petites-angle vs t"
set xlabel "t(s)"
set ylabel "angle(rad)"
set grid
plot "P7-18P-res.dat" index 0 u 1:2 w l lc rgb "orange" title "Euler cru",\
"P7-18P-res.dat" index 0 u 1:4 w l title "Predictor/corrector",\
"P7-18P-res.dat" index 1 u 1:2 w l lc rgb "#CC0605" title "Aproximació"

#----------------------------------------------------------------------------------

# Figura 2- Oscil·lacions grans. Dinàmica del pèndol.

#----------------------------------------------------------------------------------
set term png
set output "P7-18P-b-fig2.png"
set key bottom left box
set title "OSCIL·LACIONS GRANS-angle vs t"
set xlabel "temps(s)"
set ylabel "angle(rad)"
set grid
plot "P7-18P-res.dat" index 2 u 1:2 w l lc rgb "orange" title "Euler cru",\
"P7-18P-res.dat" index 2 u 1:4 w l title "Predictor/corrector"

#----------------------------------------------------------------------------------

# Figura 3- Oscil·lacions grans. Trajectòries espai fàsic.

#----------------------------------------------------------------------------------
set term png
set output "P7-18P-b-fig3.png"
set key top left box
set title "OSCIL·LACIONS GRANS-trajectòries"
set xlabel "angle(rad)"
set ylabel "velocitat(rad/s)"
set grid
plot "P7-18P-res.dat" index 2 u 2:3 w l lc rgb "orange" title "Euler cru",\
"P7-18P-res.dat" index 2 u 4:5 w l title "Predictor/corrector"

#----------------------------------------------------------------------------------

# Figura 4- Evolució energia cinètica i total.

#----------------------------------------------------------------------------------
set term png
set output "P7-18P-b-fig4.png"
set key top left box 
set title "ENERGIA CINÈTICA(K) I TOTAL (E)"
set xlabel "temps(s)"
set ylabel "energia(J)"
set grid
plot "P7-18P-res.dat" index 3 u 1:4 w l lc rgb "orange" title "E(t) euler cru ",\
"P7-18P-res.dat" index 3 u 1:7 w l title "E(t) predictor/corrector ",\
"P7-18P-res.dat" index 3 u 1:2 w l lc rgb "#EA899A" title "K(t) euler cru ",\
"P7-18P-res.dat" index 3 u 1:5 w l lc rgb "#8673A1" title "K(t) predictor/corrector "

#----------------------------------------------------------------------------------

# Figura 5- Transició. Trajectòries a l'espai fàsic.

#----------------------------------------------------------------------------------
set term png
set output "P7-18P-b-fig5.png"
set key bottom right box
set title "TRANSICIÓ - trajectòries"
set xlabel "angle(rad)"
set ylabel "velocitat (rad/s)"
set grid
plot "P7-18P-res.dat" index 4 u 2:3 w l lc rgb "#8673A1" title "alfa = 2sqrt(g/l)+0.1",\
"P7-18P-res.dat" index 4 u 4:5 w l title "alfa = 2sqrt(g/l)-0.1"

#----------------------------------------------------------------------------------

# Figura 6- Convergència del mètode. Energia total en funció del temps.

#----------------------------------------------------------------------------------
set term png
set output "P7-18P-b-fig6.png"
set key top left box
set title "CONVERGÈNCIA - E(t)"
set xlabel "temps(t)"
set ylabel "Energia total (J)"
set grid
plot "P7-18P-res.dat" index 5 u 1:2 w l lc rgb "#49678D" title "200 passos de temps",\
"P7-18P-res.dat" index 6 u 1:2 w l lc rgb "red" title "600 passos de temps",\
"P7-18P-res.dat" index 7 u 1:2 w l lc rgb "green" title "4000 passos de temps",\
"P7-18P-res.dat" index 8 u 1:2 w l lc rgb "black" title "50000 passos de temps"