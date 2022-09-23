clc; clear all; close all;

t = 0:0.1:20;

KiKp = 1e6;
dm = 1e2;
dp = 0.1;

y = (KiKp/(dm*dp))*(1 - (1/(1-dm/dp))*exp(-dm*t) - (1/(1-dp/dm))*exp(-dp*t));

figure
plot(t,y,'b', 'LineWidth', 2)
title('Producto')
xlabel('t(s)','Interpreter','latex','FontWeight','bold','FontSize',14)
ylabel('$P(\mu Mol)$','Interpreter','latex','FontWeight','bold','FontSize',14)
grid on