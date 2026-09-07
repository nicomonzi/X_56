plot [0:100] \
"output.mov" u ($0):($2) w lp title "x, solution", \
"output.echo" u ($0 + 1):($1) w p title "x, predicted", \
"output.mov" u ($0):($7*3.14159/180) w lp title "theta, solution", \
"output.echo" u ($0 + 1):($6) w p title "theta, predicted"
