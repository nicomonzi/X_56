# $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/plates/sze2004/sze31/sze31.plt,v 1.4 2010/09/11 23:56:38 masarati Exp $
#
set term post eps color blacktext "Helvetica" 22
set style line 1 lt 1 lw 2 pt 6 ps 1.5
set style line 2 lt 3 lw 2 pt 4 ps 1.5

set output 'sze31.eps'
set size 0.8,1.0
set xlabel "Tip deflections"
set ylabel "End force"
set label 1 "-Utip" at 0.5,2.0
set label 2 "Wtip" at 3.8,2.0
set key bottom right box
plot [0:7][0:4] \
"sze31_ref.dat"u 2:($1*4) title "Sze" w p ls 2, \
"sze31_16x1.dat"u 2:($1*4) title "16x1" w l ls 1, \
"sze31_ref.dat"u 3:($1*4) notitle w p ls 2, \
"sze31_16x1.dat"u 3:($1*4) notitle w l ls 1


