set term post eps enhanced color blacktext 'Helvetica' 22

set style line 1 lt 1 lw 3
set style line 2 lt 2 lw 3

set output 'asm.eps'

set key box outside

set xtic 0.01

plot \
"mov_F.txt" u 1:2 title 'x' w l ls 1, \
"mov_F.txt" u 1:3 title 'y' w l ls 2

