# $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/plates/sze2004/sze33/sze33.plt,v 1.5 2010/09/11 23:56:38 masarati Exp $
#
set term post eps color blacktext "Helvetica" 22
set style line 1 lt 1 lw 2 pt 6 ps 1.5
set style line 2 lt 3 lw 2 pt 4 ps 1.5

# for mesh in "6x30" and "10x80" ; do
#   run awk -f sze33.awk -v mesh=$mesh
#   run mbdyn -f sze33 -o sze33_$mesh
#   run ./sze33.sh sze33_$mesh
# run gnuplot sze33.plt

set output 'sze33.eps'
set size 0.8,1.0
set xlabel "Vertical deflections at point A and B (x10)"
set ylabel "Force/unit length (x0.1)"
set label 1 "WA" at 0.8,4.0
set label 2 "WB" at 1.5,4.0
set key bottom right box
plot [0:1.8][0:8] \
"sze33_ref.dat"u ($2/10):($1*8) title "Sze" w p ls 2, \
"sze33_6x30.dat"u ($2/10):($1*30) every 4 title "6x30" w p ls 1, \
"sze33_10x80.dat"u ($2/10):($1*50) title "10x80" w l ls 1, \
"sze33_ref.dat"u ($3/10):($1*8) notitle w p ls 2, \
"sze33_6x30.dat"u ($3/10):($1*30) every 4 notitle w p ls 1, \
"sze33_10x80.dat"u ($3/10):($1*50) notitle w l ls 1


