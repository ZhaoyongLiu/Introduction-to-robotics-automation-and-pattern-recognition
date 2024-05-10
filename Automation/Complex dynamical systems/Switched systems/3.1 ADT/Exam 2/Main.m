%% This MATLAB program performs a numerical simulation on Example 2.
% Link: https://zhuanlan.zhihu.com/p/539299349
% Author          Date          Version     Modification
% Zhaoyong Liu    Nov-6-2023    1.0        

%%
clc; clear; 
close all;

%% Initialization parameters
t0=0;tf=25;
t=t0:tf;   % duration
x0 = [10; 10];  % initial state

%% Switching signal
tau_a1=2;tau_a2=3;tau_a3=4;  % average dwell time
swi_seq1=tau_a1:tau_a1:tf;  % switching time sequence
swi_seq2=tau_a2:tau_a2:tf;
swi_seq3=tau_a3:tau_a3:tf;

swi_k1=Swi_signal(swi_seq1,t);
swi_k2=Swi_signal(swi_seq2,t);
swi_k3=Swi_signal(swi_seq3,t);

figure(1)
subplot(3,1,1)
stairs(t,swi_k1(2,:),'Linewidth',1.3)
axis([-inf inf 0 3])
xlabel('$k$','Interpreter','Latex')
ylabel('$\sigma$','Interpreter','Latex')
set(gca,'Linewidth',1,'Fontsize',11)

subplot(3,1,2)
stairs(t,swi_k2(2,:),'Linewidth',1.3)
axis([-inf inf 0 3])
xlabel('$k$','Interpreter','Latex')
ylabel('$\sigma$','Interpreter','Latex')
set(gca,'Linewidth',1,'Fontsize',11)

subplot(3,1,3)
stairs(t,swi_k3(2,:),'Linewidth',1.3)
axis([-inf inf 0 3])
xlabel('$k$','Interpreter','Latex')
ylabel('$\sigma$','Interpreter','Latex')
set(gca,'Linewidth',1,'Fontsize',11)

%% State trajectory
[t1,x1,x_norm1]=SLS_state([t0 tf],x0,swi_k1); % Euler method
[t2,x2,x_norm2]=SLS_state([t0 tf],x0,swi_k2);
[t3,x3,x_norm3]=SLS_state([t0 tf],x0,swi_k3);

%% Plot
% figure(2)
% plot(t,x1(1,:),'b',t,x1(2,:),'r--','Linewidth',1.3)
% grid on;
% legend('$x_1$','$x_2$','Interpreter','Latex','Fontsize',14)
% xlabel('time (s)')
% axis([-inf inf -inf inf])
% set(gca,'Linewidth',1,'Fontsize',11)

figure(3)
plot(t,x_norm1,'k--',t,x_norm2,'b-.',t,x_norm3,'r','Linewidth',1)
axis([-inf inf -1 inf])
legend('$\tau_a=2s$','$\tau_a=3s$','$\tau_a=4s$','Interpreter','Latex','Fontsize',14)
xlabel('$k$','Interpreter','Latex')
ylabel('$\|x(k)\|$','Interpreter','Latex')
set(gca,'Linewidth',1,'Fontsize',11)