%% State-dependent switched linear systems simulation 
% Reference：俞立，网络化控制系统分析与设计--切换系统处理方法[M].第47页
% Author: Liu Zhaoyong
% Date: 2023.2.13   Version: 1.0

%%
clc;clear;close all;
h=0.001;

x0=[10;10];
tspan=[0,20];
[t, x, swi_t]= SLS_Euler(@SLS_fun,tspan,x0,h); % solving differential equation

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
f=arrowPlot(x(1,:),x(2,:),'number',1);
set(f,'Color','k','LineWidth',1.2);
xlabel('$x_1$','Interpreter','Latex','Fontsize',14);
ylabel('$x_2$','Interpreter','Latex','Fontsize',14);
hold on;
x1=-40:40;y1=-0.25*x1;
x2=-20:20;y2=0.5*x2;
plot(x1,y1,'b',x2,y2,'m','linewidth',1.4);

h1=text(14,7,'$\leftarrow x_2-0.5x_1=0$');
set(h1,'Interpreter','Latex','Fontsize',13)
h2=text(-34,8.5,'$\leftarrow x_2+0.25x_1=0$');
set(h2,'Interpreter','Latex','Fontsize',13)

h3=text(-34,-2,'$\sigma=2$');
set(h3,'Interpreter','Latex','Fontsize',13)
h4=text(30,2,'$\sigma=2$');
set(h4,'Interpreter','Latex','Fontsize',13)
h5=text(0,8.2,'$\sigma=1$');
set(h5,'Interpreter','Latex','Fontsize',13)
h6=text(5,-8,'$\sigma=1$');
set(h6,'Interpreter','Latex','Fontsize',13)
box on;

%%
figure(3)
plot(t,swi_t(2,:),'linewidth',1.4,'Color','b');
xlabel('Time (s)');
legend('$\sigma(t)$','Interpreter','Latex','Fontsize',12);
axis([-inf inf 0.8 2.2]);
