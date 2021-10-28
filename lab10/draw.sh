########################################## drawing script is a given resource

set terminal post enhanced colour solid font 25

f(x, y) = 5./2. * (x**2-y)**2 + (1-x)**2

# ZALOZENIA: wyniki sa zapisane w dwoch plikach:
#   eps1.dat zawiera dwie kolumny: kolejne polozenia x_i, y_i dla epsilon = 1.e-2,
#   eps2.dat zawiera dwie kolumny: kolejne polozenia x_i, y_i dla epsilon = 1.e-3.

set xrange [-2:2]
set yrange [-2:2]
set isosample 250, 250
set table "mapa.dat"
splot f(x,y)
unset table
reset

set xrange [-2:2]
set yrange [-2:2]
set table "kontur.dat"
unset key
set contour
unset surface
set view map
set isosamples 100
set cntrparam levels 50
splot f(x,y) 
unset table

reset

set output "eps1.eps"
set xlabel "x"
set ylabel "y"
set cblabel "f(x,y)" rotate by -90
set key top center outside horizontal
set grid
set xrange [-2:2]
set yrange [-2:2]
set xtics out
set ytics out
set size ratio -1
set palette rgbformulae 33,13,10
plot "mapa.dat" with image notitle,\
 "eps1.dat" w lp lc rgb "orange" t "x_i, y_i",\
 "kontur.dat" u 1:2 w l lt -1 notitle

set output "eps2.eps"
plot "mapa.dat" with image notitle,\
 "eps2.dat" w lp lc rgb "orange" t "x_i, y_i",\
 "kontur.dat" u 1:2 w l lt -1 notitle
