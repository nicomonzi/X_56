close all
clear
clc

set(0, 'DefaultFigureWindowStyle', 'docked')

v = -3:0.001:3;
vx_Low = 1;
vx_Dmp = 0.1;

f = 0.5*(tanh((abs(v) - vx_Low)/vx_Dmp) + 1);
df = v./(2*vx_Dmp*abs(v).*cosh((abs(v) - vx_Low)/vx_Dmp).^2);

m = max(max(f), max(df));
y = [-m; m];

figure
hold on
plot(v, f, "LineWidth",2)
plot(v, df, "LineWidth",2)
plot([vx_Low, vx_Low], y, '--k')
plot(-[vx_Low, vx_Low], y, '--k')

title('$\chi$ and $\frac{\partial{\chi}}{\partial{v_x}}$, $v_{x, LOW} = 1$, $v_{x, DMP} = 0.1$', 'interpreter', 'latex')
grid on
box on
xlabel('$v_x$ $[m/s]$', 'interpreter', 'latex')
legend('$\chi$','$\frac{\partial{\chi}}{\partial{v_x}}$','interpreter','latex')
set(gca, 'TickLabelInterpreter', 'latex')
set(gca, 'FontSize', 20)

saveas(gcf, "../Plots/transition.eps", 'epsc')