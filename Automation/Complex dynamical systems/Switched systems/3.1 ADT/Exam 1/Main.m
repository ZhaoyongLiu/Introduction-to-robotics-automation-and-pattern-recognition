%% This MATLAB program performs a numerical simulation on Example 1.
% Link: https://zhuanlan.zhihu.com/p/539299349
% Author          Date          Version     Modification
% Zhaoyong Liu    Nov-3-2023    1.0        

%%
clc; clear; 
close all;

%% Initialization parameters
h = 0.001;  % step size
t0=0;tf=10;
t=t0:h:tf;   % duration
x0 = [5; 5];  % initial state

%% Switching signal
tau_a=0.5;  % average dwell time
swi_seq=0.1:tau_a:tf;  % switching time sequence
swi_t=Swi_signal(swi_seq,t);
figure(1)
plot(t,swi_t(2,:),'Color','#7F1EAC','linewidth',1.3)
axis([-inf inf 0.8 2.2])
xlabel('time (s)')
legend('$\sigma(t)$','Interpreter','Latex','Fontsize',14)
set(gca,'Linewidth',1,'Fontsize',11)

%% State trajectory
[t,x,x_norm,f]=SLS_state(@SLS_eqt,[t0 tf],x0,swi_t,h); % Euler method

%% Plot
figure(2)
plot(t,x(1,:),'b',t,x(2,:),'r--','Linewidth',1.3)
grid on;
legend('$x_1$','$x_2$','Interpreter','Latex','Fontsize',14)
xlabel('time (s)')
axis([-inf 10 -1 inf])
set(gca,'Linewidth',1,'Fontsize',11)

figure(3)
plot(t,x_norm,'k--',t,f,'r','Linewidth',1)
axis([-inf 10 -1 10])
legend('$\|x(t)\|$','$f(t)$','Interpreter','Latex','Fontsize',14)
xlabel('time (s)')
set(gca,'Linewidth',1,'Fontsize',11)