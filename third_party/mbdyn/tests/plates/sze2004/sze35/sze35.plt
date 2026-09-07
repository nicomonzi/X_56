# $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/plates/sze2004/sze35/sze35.plt,v 1.4 2010/09/13 23:40:29 masarati Exp $
#
set term post eps color blacktext "Helvetica" 22
set style line 1 lt 1 lw 2 pt 6 ps 1.5
set style line 2 lt 3 lw 2 pt 4 ps 1.5

set output 'sze35.eps'
set size 0.8,1.0
set xlabel "Displacements at points A, B and C"
set ylabel "Pulling force at point A (x10000)"
set label 1 "WA" at 1.9,2.0
set label 2 "-UB" at 2.4,1.0
set label 3 "-UC" at 3.6,1.0
set label 4 "-UB" at 4.4,1.8
set label 5 "-UC" at 3.5,3.0
set key bottom right box
plot [0:5][0:4] \
"sze35_ref.dat"u ($2):($1*4) title "Sze" w p ls 2, \
"sze35_16x24.dat"u ($2):($1/2500) every 10 title "16x24" w p ls 1, \
"sze35_24x36.dat"u ($2):($1/2500) title "24x36" w l ls 1, \
"sze35_ref.dat"u ($3):($1*4) notitle w p ls 2, \
"sze35_16x24.dat"u ($3):($1/2500) every 10 notitle w p ls 1, \
"sze35_24x36.dat"u ($3):($1/2500) notitle w l ls 1, \
"sze35_ref.dat"u ($4):($1*4) notitle w p ls 2, \
"sze35_16x24.dat"u ($4):($1/2500) every 10 notitle w p ls 1, \
"sze35_24x36.dat"u ($4):($1/2500) notitle w l ls 1


