# $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/plates/sze2004/sze36/sze36_loadinc.plt,v 1.1 2010/09/23 23:33:40 masarati Exp $
#
set term post eps color blacktext "Helvetica" 22
set style line 1 lt 1 lw 2 pt 6 ps 1.5
set style line 2 lt 3 lw 2 pt 4 ps 1.5

set output 'sze36_r.eps'
set size 0.8,1.0
set xlabel "Displacements at points A and B (x10)"
set ylabel "Pinched force, P (x10000)"
set label 1 "-WA" at 5.5,0.6
set label 2 "UB" at 1.3,0.6
set key top left width 1 box
plot [-2:9][0:1.2] \
"sze36_ref.dat"u ($2/10):($1*1.2) title "Sze" w p ls 2, \
"sze36_48x48_loadinc.dat"u ($2/10):($1/2500) title "MBDyn" w l ls 1, \
"sze36_ref.dat"u ($3/10):($1*1.2) notitle w p ls 2, \
"sze36_48x48_loadinc.dat"u ($3/10):($1/2500) notitle w l ls 1

set output 'sze36_r2.eps'
set size 0.8,1.0
set xlabel "Displacements at points A and B (x10)"
set ylabel "Pinched force, P (x10000)"
set label 1 "-WA" at 5.5,.8
set label 2 "UB" at 1.3,.8
# set key top left width 1 box
set nokey
plot [-2:10][0:5] \
"sze36_ref.dat"u ($2/10):($1*1.2) every 2 notitle w p ls 2, \
"sze36_48x48_loadinc.dat"u ($2/10):($1/2500) notitle w l ls 1, \
"sze36_ref.dat"u ($3/10):($1*1.2) every 2 notitle w p ls 2, \
"sze36_48x48_loadinc.dat"u ($3/10):($1/2500) notitle w l ls 1


