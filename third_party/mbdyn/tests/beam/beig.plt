set term post eps color blacktext "Helvetica" 22
#set term post eps color "Helvetica" 22

set output 'beig.eps'
#set key bottom right box
set key bottom left box
set xlabel "Axial compression load/critical [adim.]"
set ylabel "Bending frequency [adim.]"
#set xtic 0,.005,.02
#set ytic 0.99,.005,1.01
plot [:1][0:] \
sqrt(1-x) title "analytical" w l lt 1 lw 2, \
"beig.txt" u 1:($2/312.10) title "MBDyn, 2 el" w p lt 2 pt 5 ps 2, \
"beig.txt" u 1:($3/312.10) title "MBDyn, 4 el" w p lt 3 pt 7 ps 2
