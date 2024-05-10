%%  仿真示例
% Reference: Guang-Ren Duan, Hai-Hua Yu. LMIs in Control Systems: Analysis, Design and Applications[M]. Boca Raton: CRC Press, 2013, Page 201.
% Author          Date           Version     Modification
% Zhaoyong Liu    May-9-2024     1.0         change ode4 with ode45

%%
clc; clear;
close all;

%% Initialization parameters
t0 = 0;  tf = 5;
x0 = (2:-0.1:0.1)'; % initial state
global K;
load('K.mat','K')

%% State trajectories of the system 
[t,x]=ode45(@sys_func,[t0 tf],x0);  

%% Plot the figure
figure(1)
plot(t,x,'Linewidth',1.3)
xlabel('time (s)')
set(gca,'Linewidth',1,'Fontsize',11)


%% System function
function dotx = sys_func(t,x)
% input: time, state
% output: derivative of the state

n = 20;
e1 = 20*ones(n,1);
e2 = (20: -1: 1)';
Asp = spdiags([e1 e2],-1:0,n,n);
A = full(Asp);
B = [1 zeros(1,n-1)]';

global K;

u = K*x; % control input
dotx = A*x+B*u;

end
