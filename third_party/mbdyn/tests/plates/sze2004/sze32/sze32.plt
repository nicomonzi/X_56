# $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/plates/sze2004/sze32/sze32.plt,v 1.3 2010/09/11 23:56:38 masarati Exp $
#
set term post eps color blacktext "Helvetica" 22
set style line 1 lt 1 lw 2 pt 6 ps 1.5
set style line 2 lt 3 lw 2 pt 4 ps 1.5

set output 'sze32.eps'
set size 0.8,1.0
set xlabel "Tip deflections (x10)"
set ylabel "End moment (x50/3)"
set label 1 "-Utip" at 1.5,0.6
set label 2 "Wtip" at 0.2,0.6
set key bottom right box
plot [0:2][0:1] \
"sze32_ref.dat"u ($2/10):($1) title "Sze" w p ls 2, \
"sze32_16x1.dat"u ($2/10):($1) title "16x1" w l ls 1, \
"sze32_ref.dat"u ($4/10):($1) notitle w p ls 2, \
"sze32_16x1.dat"u ($3/10):($1) notitle w l ls 1


