set term post eps color blacktext "Helvetica" 22
#set term post eps color "Helvetica" 22

set output 'beig_tf.eps'
#set key bottom right box
set key bottom left width -6 box
set xlabel "Angular frequency [radian/s]"
set ylabel "Displacement [m/N]"
#set xtic 10,2.5,1000
#set ytic 0.99,.005,1.01
set logscale xy
plot [10:4000] \
"beig_mbdyn.txt" title "MBDyn" w p lt 1 pt 7 ps .5, \
"beig_rom.txt" title "ROM - 1 mode" w l lt 2 lw 3, \
"beig_rom2.txt" title "ROM - 2 modes" w l lt 3 lw 3, \
"beig_rom3.txt" title "ROM - 1 mode + residualization" w l lt 4 lw 3

