%% Trajectories of output and estimated output
% Reference: Guang-Ren Duan, Hai-Hua Yu. LMIs in Control Systems: Analysis, Design and Applications[M]. Boca Raton: CRC Press, 2013, Page 306.
% Author          Date           Version     Modification
% Zhaoyong Liu    Jun-20-2024     1.0   

%%
clc; clear;
close all;

%% Initialization parameters
t0 = 0;  tf = 6;
h = 0.001;  % sampling period
x0 = [3 2 1]'; % initial state
dimn = 3;
zeta0 = zeros(dimn,1);
extendedState0 = [x0; zeta0];

global Af Bf Cf Df;
load('Af.mat','Af');
load('Bf.mat','Bf');
load('Cf.mat','Cf');
load('Df.mat','Df');

%% State trajectories of system and filter
[t,extendedState,zHatt] = euler(@sys_flt_func,[t0 tf],extendedState0,h); 

x = extendedState(1:dimn,:);
zeta = extendedState(dimn+1:end,:);

L =[-1.1604 1.7916 -1.4033;
    -2.0351 1.5369 -0.3777];
z = L*x; % interested output

%% Plot the figures
figure(1)
subplot(2,1,1)
plot(t,z(1,:),'r',t,zHatt(1,:),'b--','Linewidth',1.3)
legend('$z_1$','$\hat{z}_1$','Interpreter','Latex','Fontsize',14)
xlabel('time (s)')
axis([-inf inf -inf inf])
set(gca,'Linewidth',1,'Fontsize',11)
subplot(2,1,2)
plot(t,z(2,:),'r',t,zHatt(2,:),'b--','Linewidth',1.3)
legend('$z_2$','$\hat{z}_2$','Interpreter','Latex','Fontsize',14)
xlabel('time (s)')
axis([-inf inf -inf inf])
set(gca,'Linewidth',1,'Fontsize',11)

figure(2)
subplot(2,1,1)
plot(t,z(1,:)-zHatt(1,:),'k','Linewidth',1.3)
legend('$\tilde{z}_1$','Interpreter','Latex','Fontsize',14)
xlabel('time (s)')
axis([-inf inf -inf inf])
set(gca,'Linewidth',1,'Fontsize',11)
subplot(2,1,2)
plot(t,z(2,:)-zHatt(2,:),'k','Linewidth',1.3)
legend('$\tilde{z}_2$','Interpreter','Latex','Fontsize',14)
xlabel('time (s)')
axis([-inf inf -inf inf])
set(gca,'Linewidth',1,'Fontsize',11)

%% System and filter functions
function [dotExtendedState,zHat] = sys_flt_func(t,extendedState)
% input arguments: time, state, filter state
% output arguments: derivatives of state and filter state, estimated output zHat(t)

A = [-4.4697 9.3324  8.7645;
     -9.6408 2.0236  2.2076;
     -6.6937 1.3595 -3.2241];
B = [0.2309;
    -0.4584;
    -0.8390];
C = [ 0.7813 -4.9301 -4.3341;
      4.9421 -3.4866 -4.1502];
D = [0.0074;
     -0.0088];

global Af;
global Bf;
global Cf;
global Df;

n = 3;
x = extendedState(1:n);
zeta = extendedState(n+1:end);
w = wgn(1,1,10); 
dotx = A*x+B*w; 
y = C*x+D*w; 
dotzeta = Af*zeta+Bf*y; 
zHat = Cf*zeta+Df*y;
dotExtendedState = [dotx; dotzeta];

end

%% Forward euler method
function [t,extendedState,zHatt]=euler(ufunc,tspan,extendedState0,h)
    if nargin<4
        h=0.01;
    end
    if length(tspan)==2
        t0=tspan(1);
        tn=tspan(2);
    else
        error(message('MATLAB: Wrong Dimension Of Tspan'));
    end
    n=floor((tn-t0)/h); 
    t=zeros(1,n+1);
    t(1)=t0; 
    extendedState(:,1)=extendedState0;
    zHatt = zeros(2,n+1);
    for i=1:n
        t(i+1)=t(i)+h;
        [dotExtendedState, zHat] = ufunc(t(i),extendedState(:,i));
        extendedState(:,i+1)=extendedState(:,i)+h*dotExtendedState;
        zHatt(:,i) = zHat;
    end
end

