set term post eps color blacktext "Arial" 22
set grid

set style line 1 lt 1 lw 4 pt 5 ps 2
set style line 2 lt 2 lw 4 pt 7 ps 2
set style line 3 lt 3 lw 4 pt 9 ps 2

set output 'torque-1.eps'
set nokey
set xlabel "Time, s"
set ylabel "Torque, adim"
p [0:10][-1.1:1.1] \
"/tmp/modal.frc" u (1e-3*$0):(-$5/4.85) notitle w l ls 1

set output 'tip-y.eps'
set key bottom right box
set xlabel "Time, s"
set ylabel "Displacement, m"
p [0:10][-.01:.09] \
"<awk '$1==2000' /tmp/modal.mov" u (1e-3*$0):($3) title "modal" w l ls 1, \
"<awk '$1==1010' /tmp/beam10.mov" u (1e-3*$0):($3) title "beam, 10 nodes" w l ls 2, \
"<awk '$1==1020' /tmp/beam20.mov" u (1e-3*$0):($3) title "beam, 20 nodes" w l ls 3

set output 'root-a.eps'
set key bottom left box
set xlabel "Time, s"
set ylabel "rotation, deg"
p [0:10][:] \
"<awk '$1==1000' /tmp/modal.mov" u (1e-3*$0):($7) title "modal" w l ls 1, \
"<awk '$1==1000' /tmp/beam10.mov" u (1e-3*$0):($7) title "beam, 10 nodes" w l ls 2, \
"<awk '$1==1000' /tmp/beam20.mov" u (1e-3*$0):($7) title "beam, 20 nodes" w l ls 3

set output 'axial.eps'
set nokey
set xlabel "Time, s"
set ylabel "Axial force, N"
set ytic -300,100,300
set label 1 "buckling load" at 0.2,-115
p [0:10][-250:250] \
"<awk '$1==1' /tmp/beam10.act" u (1e-3*$0):($2) notitle w l ls 1, \
-87 notitle w l lt -1 lw 4


