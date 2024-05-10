%% Solving Difference Equations

function [t,x,x_norm]=SLS_state(tspan,x0,swi_k)
% Input arguments: time span, initial state, switching signal
% Output arguments: time, state, norm of state x,

t0=tspan(1);
tf=tspan(2);
n=tf-t0; 

t=zeros(1,n+1);t(1)=t0;  
x=zeros(2,n+1);x(:,1)=x0;  % system state
x_norm=zeros(1,n+1);x_norm(:,1)=norm(x0);   % norm of state x

A1=[-0.8 0.2;-0.5 -0.3]; % system matrices
A2=[-0.4 0.5;-0.2 -0.5];

for i=1:n
    t(i+1)=t(i)+1;
    switch swi_k(2,i)   
       case 1
           A=A1;
       case 2
           A=A2;
    end
    x(:,i+1)=A*x(:,i);  % update the state
    x_norm(:,i+1)=norm(x(:,i+1));
end