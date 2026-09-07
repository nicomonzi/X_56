# $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/plates/chimakurthi2009/case4/case4.plt,v 1.3 2011/10/19 09:57:31 masarati Exp $
#
set term post eps color blacktext "Helvetica" 22
set style line 1 lt 1 lw 2
set style line 2 lt 3 lw 2 pt 4 ps 1.5
set style line 3 lt 4 lw 2 pt 6 ps 1.5
set style line 4 lt 2 lw 2

set output 'case4_05.eps'
set xlabel "Time, t/T"
set ylabel "Tip displacement, z/L"
set key bottom right box
plot [][-.2:.7] \
"y_05_32000"u 1:2 every 10 title "32x9" w l ls 1, \
"z_05_i.dat"u 1:2 title "MSC.Marc" w p ls 2


set output 'case4_10.eps'
set xlabel "Time, t/T"
set ylabel "Tip displacement, z/L"
set key bottom right box
plot [][-.2:.7] \
"yy_10_16000"u 1:2 every 10 title "16x6" w l ls 1, \
"z_10_i.dat"u 1:2 title "UM/NLAMS" w p ls 2, \
"z_10_marc_i.dat"u 1:2 title "MSC.Marc" w p ls 3


set output 'case4_30.eps'
set xlabel "Time, t/T"
set ylabel "Tip displacement, z/L"
set key bottom right box
plot [:4][-1.2:1.2] \
"y_30_32000"u 1:2 every 50 title "32x9" w l ls 1, \
"z_30_i.dat"u 1:2 title "MSC.Marc" w lp ls 2


