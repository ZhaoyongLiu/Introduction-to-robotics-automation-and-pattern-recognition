%% Simulation of state-dependent switched linear systems
% Reference: 俞立, 网络化控制系统分析与设计--切换系统处理方法[M] 第46页
% Author: Liu Zhaoyong
% Date: 2023.1.13   Version: 1.0

%%
clc;clear;close all;
h=0.001; % step size

x0=[10;10];
tspan=[0,30];
[t, x]= SLS_Euler(@SLS_fun,tspan,x0,h); % solving differential equation

%% plot
figure(1)
hold on;
plot(t,x(1,:),'k-','Linewidth',1.5);
plot(t,x(2,:),'Linestyle','--','Color','#77AC30','Linewidth',1.5);
xlabel('Time(s)')
legend('$x_1$','$x_2$','Interpreter','latex','Fontsize',13);
box on;
%%
figure(2)
hold on;
f=arrowPlot(x(1,:),x(2,:),'number',1);
set(f,'LineWidth',1.2);
xlabel('$x_1$','Interpreter','Latex','Fontsize',14);
ylabel('$x_2$','Interpreter','Latex','Fontsize',14);
box on;