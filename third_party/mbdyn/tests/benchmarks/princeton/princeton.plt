set term post eps color enhanced blacktext 'Helvetica' 22

set style line 1 lt 1 lw 4 pt 5 ps 2
set style line 2 lt 2 lw 4 pt 7 ps 2
set style line 3 lt 3 lw 4 pt 9 ps 2

set output 'princeton_u2.eps'

set xlabel 'Loading angle, {/Symbol q} [DEG]'
set ylabel 'Displacement, u_2 [m]'
set xtic 0,15,90
set key top left box width 1

plot [0:90][0:.26] \
"princeton_p1.dat" u 1:2 title 'P1' w lp ls 1, \
"princeton_p2.dat" u 1:2 title 'P2' w lp ls 2, \
"princeton_p3.dat" u 1:2 title 'P3' w lp ls 3



set output 'princeton_u3.eps'

set xlabel 'Loading angle, {/Symbol q} [DEG]'
set ylabel 'Displacement, u_3 [m]'
set xtic 0,15,90
set key top right box width 1

plot [0:90][0:.021] \
"princeton_p1.dat" u 1:3 title 'P1' w lp ls 1, \
"princeton_p2.dat" u 1:3 title 'P2' w lp ls 2, \
"princeton_p3.dat" u 1:3 title 'P3' w lp ls 3



set output 'princeton_phi.eps'

set xlabel 'Loading angle, {/Symbol q} [DEG]'
set ylabel 'Twist angle, {/Symbol f} [RAD]'
set xtic 0,15,90
set key top right box width 1

plot [0:90][0:.07] \
"princeton_p1.dat" u 1:4 title 'P1' w lp ls 1, \
"princeton_p2.dat" u 1:4 title 'P2' w lp ls 2, \
"princeton_p3.dat" u 1:4 title 'P3' w lp ls 3

