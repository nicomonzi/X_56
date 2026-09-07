set term post eps color blacktext 'Helvetica' 22

set style line 1 lt -1 lw 1
set style line 2 lt 2 lw 2
set style line 3 lt 3 lw 2
set style line 4 lt 1 lw 2

set output 'o0_10.eps'

set key top right box width -4
set xlabel "Time, s"
set ylabel "Tip displacement, m"

plot [0:10] \
"o0_10_yoo.dat" title "Yoo & Chung, 2001" w l ls 1, \
"o0_10_jinyang.dat" title "Jinyang & Jiazhen, 2005" w l ls 2, \
"o0_10_singh.dat" u 1:(-$2) title "Singh & Chopra, 2006" w l ls 3, \
"o0_10_32x64_present.dat" title "present, 32x64" w l ls 4

set output 'o0_10_mesh_refinement.eps'

set key top right box width -2
set xlabel "Time, s"
set ylabel "Tip displacement, m"

plot [0:10] \
"o0_10_8x16_present.dat" title "present, 8x16" w l ls 2, \
"o0_10_16x32_present.dat" title "present, 16x32" w l ls 3, \
"o0_10_32x64_present.dat" title "present, 32x64" w l ls 4

