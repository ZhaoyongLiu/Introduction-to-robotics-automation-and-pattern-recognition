%% matlab delayed-diff system simulation
clear; close all; clc;
lags=0.4;
t_vec=(0:1e-4:20)';
fig1=figure(1);fig1.Color=[1,1,1];
for ii=1:1:1
    sol1=dde23(@dde_func,lags,@x_history,t_vec);
    x_vec = deval(sol1,t_vec);
    x1_trajectory=x_vec(1,:);
    x2_trajectory=x_vec(2,:);
%     plot(t_vec,x1_trajectory,'LineStyle','-','LineWidth',1.5,'Color','r');hold on;
%     plot(t_vec,x2_trajectory,'LineStyle','-','LineWidth',1.5,'Color','b');
%     legend('x1','x2');xlabel('time');ylabel('x(t)');
    t_vec_minus=(-lags:1e-4:0-1e-4)';
    t_vec_whl=[t_vec_minus;t_vec];
    x_vec_minus=[sin(t_vec_minus)';sin(2*t_vec_minus)'];
    x_vec_whl=[x_vec_minus x_vec];
    x1_trajectory_whl=x_vec_whl(1,:);
    x2_trajectory_whl=x_vec_whl(2,:);
    plot(t_vec_whl,x1_trajectory_whl,'LineStyle','-','LineWidth',1.6,'Color','r');hold on;
    plot(t_vec_whl,x2_trajectory_whl,'LineStyle','-','LineWidth',1.6,'Color','b');
    axis([-lags inf -inf inf]);
    h1=legend('$x_1(t)$','$x_2(t)$');set(h1,'Interpreter','latex','FontSize',13);
    xlabel('Time (s)');
    h1=ylabel('$x(t)$');set(h1,'Interpreter','latex','FontSize',13);
end

function xdot=dde_func(t,x,x_delayed)
x1=x(1);
x2=x(2);
xdot=zeros(2,1);
x1_delayed=x_delayed(1);
x2_delayed=x_delayed(2);
   A=[-2 0;0 -0.9];
   Ad=0.899*[-1 0;-1 -1];
   xdot=A*[x1;x2]+Ad*[x1_delayed;x2_delayed];
end

function x=x_history(t)
%    x=0.1*ones(2,1)+rand(2,1)*0.3;
   x=[sin(t);sin(2*t)];
%    x=0.1*ones(2,1);
end