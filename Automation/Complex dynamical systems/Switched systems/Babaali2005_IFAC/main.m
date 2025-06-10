%% State and switching signal estimation for discrete-time linear switched systems
% Reference: Babaali, et al. ASYMPTOTIC OBSERVERS FOR DISCRETE-TIME SWITCHED LINEAR
% SYSTEMS[C]. IFAC, 2005.

% Author         Date         Version    Modification
% Zhaoyong Liu   Jun-9-2025   1.0

%%
clc; clear; close all;

A1 = [1 0; 0 0.5]; C1 = [1 0];
A2 = [2 0; 0 0.5]; C2 = [2 0];

L1 = [0.5; 0]; L2 = [1.5; 0];

A = {A1,A2}; C = {C1,C2};
L = {L1,L2};

matrices = {A,C,L};

x0 = [1; 1];
k0 = 0; kf = 10;
swiSigHat0 = 1;
xhat0 = [0; 2];

initialValues = {x0, swiSigHat0, xhat0};

%% Switching signal 
swiSig = zeros(2, kf-k0+1);
for i = k0:1:kf
    swiSig(1,i+1) = i;
    if mod(i,2) == 1
        swiSig(2,i+1) = 1;
    elseif mod(i,2) == 0
        swiSig(2,i+1) = 2;
    end
end

%% Simulation
[t,x, xhat, swiSigHat] = JointObsvDSLS(@DSLS,@JointObsv,[k0 kf],swiSig,initialValues,matrices);

%% Plot the figures
figure(1)

plot(t,x(1,:),'bsquare',t,xhat(1,:),'r*','MarkerSize',10,'LineWidth',1.2)
legend('$x_1$','$\hat{x}_1$','Interpreter','latex','Fontsize',16)
xlabel('Time')

figure(2)

plot(t,x(2,:),'square','Color',"#D95319",'MarkerSize',10,'LineWidth',1.2)
hold on;
plot(t,xhat(2,:),'*','COlor',"#77AC30",'MarkerSize',10,'LineWidth',1.2)
legend('$x_2$','$\hat{x}_2$','Interpreter','latex','Fontsize',16)
xlabel('Time')

figure(3)
plot(t,swiSig(2,:),'s','Color',"#7E2F8E",'MarkerSize',10,'LineWidth',1.2)
hold on;
plot(t,swiSigHat(2,:),"x",'Color','r','MarkerSize',10,'LineWidth',1.2)
legend('$\sigma$','$\hat{\sigma}$','Interpreter','latex','Fontsize',16)
xlabel('Time')
axis([-inf inf 0.5 2.5])
