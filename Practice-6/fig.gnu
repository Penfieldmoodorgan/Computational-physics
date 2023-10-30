#--------------------------------------------------------------------------------------------

# Figura 1- Convergència del càlcul de I1 dibuixant l'error estimat comparant-lo amb el real

#---------------------------------------------------------------------------------------------
set term png
set output "P6-18P-fig1.png"
set key top right
set xlabel "n"
set ylabel "I1"
set grid
set format x '%.e'
f=(122./27.)-(21./4.)*(pi**2)
plot "P6-18P-res.dat" index 0 u 1:2 w l title "I1 MC",\
"P6-18P-res.dat" index 0 u 1:2:3 w yerrorbars lc rgb "#C51D34" notitle,\
f w l title "exacte"

#-------------------------------------------------------------------------------------

# Figura 2- Convergència del càlcul de I2 

#----------------------------------------------------------------------------------------
set term png
set output "P6-18P-fig2.png"
set key top right
set xlabel "n"
set ylabel "I2"
set grid
set format x '%.e'
plot "P6-18P-res.dat" index 2 u 1:2 w l title "I2 MC",\
"P6-18P-res.dat" index 2 u 1:2:3 w yerrorbars notitle

#-------------------------------------------------------------------------------------

# Figura 3- Convergència del càlcul de I3 

#----------------------------------------------------------------------------------------
set term png
set output "P6-18P-fig3.png"
set key top right
set xlabel "n"
set ylabel "I3"
set grid
set format x '%.e'
plot "P6-18P-res.dat" index 3 u 1:2 w l title "I3 MC",\
"P6-18P-res.dat" index 3 u 1:2:3 w yerrorbars notitle
