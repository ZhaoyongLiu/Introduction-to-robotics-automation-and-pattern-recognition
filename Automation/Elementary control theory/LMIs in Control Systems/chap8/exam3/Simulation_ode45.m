%% Trajectories of state and estimated state (ode45)
% Reference: Guang-Ren Duan, Hai-Hua Yu. LMIs in Control Systems: Analysis, Design and Applications[M]. Boca Raton: CRC Press, 2013, Page 298.
% Author          Date           Version     Modification
% Zhaoyong Liu    Jun-17-2024     1.0   

%%
clc; clear;
close all;
tic;

%% Initialization parameters
t0 = 0;  tf = 10;
x0 = [3 2 1]'; % initial state
dimn = 3;
xHat0 = zeros(dimn,1);
extendedState0 = [x0; xHat0];

global L;
load('L.mat','L');


%% State trajectories of system and observer
[t,extendedState] = ode45(@sys_obv_func,[t0 tf],extendedState0);  
extendedState = extendedState';
x = extendedState(1:dimn,:);
xHat = extendedState(dimn+1:end,:);

C2 = [0.9607 1.5600 2.8558;
     -2.4371 1.3634 0.0095];
z = C2*x; % output
zHat = C2*xHat;

%% Plot the figures
figure(1)
plot(t,z(1,:)-zHat(1,:),'k','Linewidth',1.3)
legend('$\tilde{z}_1$','Interpreter','Latex','Fontsize',14)
xlabel('time (s)')
axis([-inf inf -inf inf])
set(gca,'Linewidth',1,'Fontsize',11)


figure(2)
plot(t,z(2,:)-zHat(2,:),'k','Linewidth',1.3)
legend('$\tilde{z}_2$','Interpreter','Latex','Fontsize',14)
xlabel('time (s)')
axis([-inf inf -inf inf])
set(gca,'Linewidth',1,'Fontsize',11)

figure(3)
plot(t,x(1,:),t,x(2,:),'--',t,x(3,:),'-.','Linewidth',1.3)
legend('$x_1$','$x_2$','$x_3$','Interpreter','Latex','Fontsize',14)
xlabel('time (s)')
axis([-inf inf -inf inf])
set(gca,'Linewidth',1,'Fontsize',11)
toc;

%% System and observer functions
function dotExtendedState = sys_obv_func(t,extendedState)
% input: time, state, estimated state
% output: derivatives of state and estimated state

A = [2.8982 -3.1606  0.6816;
     6.3595 -4.2055  4.5423;
     3.2046 -3.1761 -3.8142];
B1 = [2.0310  1.2164;
      0.4084  0.2794;
     -0.7775 -0.3307];
B2 = [-0.0785; -0.0853; 0.0986];
C1 = [-0.8778 -4.9442 -4.5084;
       4.0161 -2.0259  1.9318];
C2 = [0.9607 1.5600 2.8558;
     -2.4371 1.3634 0.0095];
D1 = [0.6004  0.2107;
      1.9320 -0.3997];
D2 = [0.0330; -0.0414];
global L;

n = 3;
x = extendedState(1:n);
xHat = extendedState(n+1:end);
u = [2*sin(5*t); cos(2*t)]; % control input
w = wgn(1,1,10); 
dotx = A*x+B1*u+B2*w; 
y = C1*x+D1*u+D2*w; 
dotxHat = (A-L*C1)*xHat+L*y+(B1-L*D1)*u; 
dotExtendedState = [dotx; dotxHat];

end
