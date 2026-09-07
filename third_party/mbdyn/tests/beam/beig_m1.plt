set term post eps color blacktext "Helvetica" 22
#set term post eps color "Helvetica" 22

set output 'beig_m1.eps'
#set key bottom right box
#set key bottom left box
set key .57,.15 box
set xlabel "Beam span [adim.]"
set ylabel "Beam deflection [adim.]"
#set xtic 0,.005,.02
#set ytic 0.99,.005,1.01
plot [:1][0:] \
sin(3.14159*x) title "analytical" w l lt 1 lw 2, \
"beig_m1.txt" u 1:3 title "MBDyn, 2 el" w p lt 2 pt 5 ps 2, \
"beig_m1.txt" u 1:2 title "MBDyn, 4 el" w p lt 3 pt 7 ps 2
