%% 
% Author         Date         Version    Modification
% Zhaoyong Liu   Jun-10-2025  1.0

%%
clc; clear; close all;

h = 0.001;   % step size
t0 = 0; tf = 20;
t = t0:h:tf;

%% Control input
SwPer=5;
u=rem(floor(t/SwPer),2)-0.5;

%% System Matrices
A={[0 1;-1 0];[0 1;-1.1 0];[0 1;-1.2 0];[0 1;-1.25 0]};
B=[0; 1];

%% Numerical simulation
x0 = [0.01; 0];
x = zeros(2,tf/h);
x(:,1) = x0;
tic
for k = 1:tf/h
    x(:,k+1) = x(:,k) + h*(A{sigma(t(k))}*x(:,k) + B*u(k));
end
toc

%% Phase Portrait
figure(1)
plot(x(1,:),x(2,:),'LineWidth',2)
xlabel({'$x_1(t)$'},'Interpreter','latex') 
ylabel({'$x_2(t)$'},'Interpreter','latex') 
set(gca,'fontsize',18)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);

%% Time response
figure(2)
plot(t,x(1,:),t,x(2,:),'LineWidth',2)
legend('$x_1(t)$','$x_2(t)$','Interpreter','latex')
set(gca,'fontsize',18)

%% Switching signal
function i = sigma(t)
if (t<5)
    i=1;
elseif (5<=t&&t<10)
    i=2;
elseif (10<=t&&t<15)
    i=3;
elseif (15<=t)
    i=4;
end
end

