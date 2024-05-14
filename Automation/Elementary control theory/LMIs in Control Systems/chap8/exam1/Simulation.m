%% Trajectories of state and estimated state
% Reference: Guang-Ren Duan, Hai-Hua Yu. LMIs in Control Systems: Analysis, Design and Applications[M]. Boca Raton: CRC Press, 2013, Page 287.
% Author          Date           Version     Modification
% Zhaoyong Liu    May-9-2024     1.0   

%%
clc; clear; 
close all;

%% Initialization parameters
t0 = 0;  tf = 20;
x0 = [3 2 1]'; % initial state
dimn = 3;
xhat0 = zeros(dimn,1);
extendedState0 = [x0; xhat0];

global L;
load('L.mat','L');

%% State trajectories of system and observer
[t,extendedState]=ode45(@sys_obv_func,[t0 tf],extendedState0);  
extendedState = extendedState';
x = extendedState(1:dimn,:);
xhat = extendedState(dimn+1:end,:);

%% Plot the figures
figure(1)
subplot(3,1,1)
plot(t,x(1,:)-xhat(1,:),'k','Linewidth',1.3)
legend('$e_1$','Interpreter','Latex','Fontsize',14)
xlabel('time (s)')
axis([-inf inf -inf inf])
set(gca,'Linewidth',1,'Fontsize',11)

subplot(3,1,2)
plot(t,x(2,:)-xhat(2,:),'k','Linewidth',1.3)
legend('$e_2$','Interpreter','Latex','Fontsize',14)
xlabel('time (s)')
axis([-inf inf -inf inf])
set(gca,'Linewidth',1,'Fontsize',11)

subplot(3,1,3)
plot(t,x(3,:)-xhat(3,:),'k','Linewidth',1.3)
legend('$e_3$','Interpreter','Latex','Fontsize',14)
xlabel('time (s)')
axis([-inf inf -inf inf])
set(gca,'Linewidth',1,'Fontsize',11)

figure(2)
subplot(3,1,1)
plot(t,x(1,:),'k',t,xhat(1,:),'r--','Linewidth',1.3)
legend('$x_1$','$\hat{x}_1$','Interpreter','Latex','Fontsize',14)
xlabel('time (s)')
axis([-inf inf -inf inf])
set(gca,'Linewidth',1,'Fontsize',11)

subplot(3,1,2)
plot(t,x(2,:),'k',t,xhat(2,:),'r--','Linewidth',1.3)
legend('$x_2$','$\hat{x}_2$','Interpreter','Latex','Fontsize',14)
xlabel('time (s)')
axis([-inf inf -inf inf])
set(gca,'Linewidth',1,'Fontsize',11)

subplot(3,1,3)
plot(t,x(3,:),'k',t,xhat(3,:),'r--','Linewidth',1.3)
legend('$x_3$','$\hat{x}_3$','Interpreter','Latex','Fontsize',14)
xlabel('time (s)')
axis([-inf inf -inf inf])
set(gca,'Linewidth',1,'Fontsize',11)

%% System and observer functions
function dotExtendedState = sys_obv_func(t,extendedState)
% input: time, state, estimated state
% output: derivatives of state and estimated state

A = [-0.5  0  0;
        0 -2 10;
        0  1 -2 ];
B = [ 1 0;
     -2 2;
      0 1];
C = [1 0 0;
     0 0 1];
global L;

n = 3;
x = extendedState(1:n);
xhat = extendedState(n+1:end);
u = [2*sin(5*t)+0.1*t; cos(2*t)]; % control input
dotx = A*x+B*u; 
y = C*x;
dotxhat = A*xhat+B*u+L*(y-C*xhat); 
dotExtendedState = [dotx; dotxhat];

end
