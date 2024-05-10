%% Forward Euler Method for Solving Differential Equations

function [t,x,x_norm,f]=SLS_state(SLS_eqt,tspan,x0,swi_t,h)
% Input arguments: function name, time span, initial state, switching signal,step size
% Output arguments: time, state, norm of state x, f(t)

t0=tspan(1);
tf=tspan(2);
n=floor((tf-t0)/h); 

t=zeros(1,n+1);t(1)=t0;  
x=zeros(2,n+1);x(:,1)=x0;  % system state
x_norm=zeros(1,n+1);x_norm(:,1)=norm(x0);   % norm of state x

% estimatiing function f satisfying ||x(t)||<= f(t) = sqrt(epsilon2/epsilon1)*exp(-alpha(t-t0))*||x0||
P1=[5.7363 0;0 2.8702]; P2=[2.8702 0;0 5.7363];
lambda_min1=min(eig(P1)); lambda_min2=min(eig(P2));
lambda_max1=max(eig(P1)); lambda_max2=max(eig(P2));

epsilon1=min(lambda_min1,lambda_min2);
epsilon2=max(lambda_max1,lambda_max2);

lambda=3.999; mu=1.9986; tau_a=0.5;
alpha=1/2*(lambda-log(mu)/tau_a);

f=zeros(1,n+1); f(:,1)=sqrt(epsilon2/epsilon1)*norm(x0);

for i=1:n
    t(i+1)=t(i)+h;
    dx=SLS_eqt(swi_t(2,i),x(:,i));
    x(:,i+1)=x(:,i)+h*dx;  % update the state
    x_norm(:,i+1)=norm(x(:,i+1));
    f(:,i+1)=sqrt(epsilon2/epsilon1)*exp(-alpha*t(i+1))*norm(x0);
end