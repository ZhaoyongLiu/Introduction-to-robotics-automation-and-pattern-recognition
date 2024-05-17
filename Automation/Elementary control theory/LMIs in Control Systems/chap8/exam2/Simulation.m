%% Trajectories of state and estimated state
% Reference: Guang-Ren Duan, Hai-Hua Yu. LMIs in Control Systems: Analysis, Design and Applications[M]. Boca Raton: CRC Press, 2013, Page 292.
% Author          Date           Version     Modification
% Zhaoyong Liu    May-17-2024     1.0   

%%
clc; clear; 
close all;

%% Initialization parameters
t0 = 0;  tf = 20;
x0 = [3 2 1]'; % initial state
dimn = 3;
z0 = zeros(1,1);
extendedState0 = [x0; z0];

global F G H M N;
load('F.mat','F');
load('G.mat','G');
load('H.mat','H');
load('M.mat','M');
load('N.mat','N');

%% State trajectories of system and observer
[t,extendedState]=ode45(@sys_obv_func,[t0 tf],extendedState0);  
extendedState = extendedState';
x = extendedState(1:dimn,:);
z = extendedState(dimn+1:end,:);

C = [1 0 0;
     0 0 1];
y = C*x; % output
xi2Hat = M*z+N*y;

%% Plot the figures
figure(1)
plot(t,x(2,:)-xi2Hat(:,:),'k','Linewidth',1.3)
legend('$e_2$','Interpreter','Latex','Fontsize',14)
xlabel('time (s)')
axis([-inf inf -0.5 inf])
set(gca,'Linewidth',1,'Fontsize',11)


figure(2)
plot(t,x(2,:),'k',t,xi2Hat(:,:),'r--','Linewidth',1.3)
legend('$x_2$','$\hat{x}_2$','Interpreter','Latex','Fontsize',14)
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
global F G H;

n = 3;
x = extendedState(1:n);
z = extendedState(n+1:end);
u = [2*sin(5*t)+0.1*t; cos(2*t)]; % control input
dotx = A*x+B*u; 
y = C*x;
dotz = F*z+G*y+H*u; 
dotExtendedState = [dotx; dotz];

end
