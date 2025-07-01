%%
clc; clear; close all;
t0=1; tn=2;
tspan=[t0,tn];
y0=2/5;
h=0.1;

[t,y]=anal(tspan,y0,h);
[t1,y1]=euler(@ode,tspan,y0,h);
[t2,y2]=backeuler(tspan,y0,h);
[t3,y3]=heun(@ode,tspan,y0,h);
[t4,y4]=ralston(@ode,tspan,y0,h);
[t5,y5]=runge_kutta_4(@ode,tspan,y0,h);

% plot
figure(1)
hold on;
plot(t,y,'square-','Color','k','Linewidth',1.2);
plot(t1,y1,'--','Linewidth',1.2);
plot(t2,y2,'-.','Linewidth',1.2);
plot(t3,y3,'o-','Linewidth',1.2);
plot(t4,y4,'*-','Linewidth',1.2);
plot(t5,y5,'diamond-','Linewidth',1.2);
legend('Analytical solution','Forward euler','Backward euler','Heun','Ralston','Runge-Kutta 4','Fontsize',12);
hold off;
xlabel('Time/s')
axis([-inf tn -inf inf])

%% Dynamic equation
function dy=ode(t,y)
    dy=t^3-y/t;
end

%% Analytical solution
function [t,y]=anal(tspan,y0,h)
    if length(tspan)==2
        t0=tspan(1);
        tn=tspan(2);
    else
        error('MATLAB: Wrong Dimension Of Tspan');
    end
    n=floor((tn-t0)/h); 
    t=zeros(1,n+1);
    t(1)=t0; 
    y(:,1)=y0;
    for i=1:n
        t(i+1)=t(i)+h;
        y(:,i+1)=1/5*t(i+1)^4+1/(5*t(i+1));
    end
end

%% Forward euler method
function [t,y]=euler(ufunc,tspan,y0,h)
    if nargin<4
        h=0.01;
    end
    if length(tspan)==2
        t0=tspan(1);
        tn=tspan(2);
    else
        error('MATLAB: Wrong Dimension Of Tspan');
    end
    n=floor((tn-t0)/h); 
    t=zeros(1,n+1);
    t(1)=t0; 
    y(:,1)=y0;
    for i=1:n
        t(i+1)=t(i)+h;
        y(:,i+1)=y(:,i)+h*ufunc(t(i),y(:,i));
    end
end

%% Backward euler method
function [t,y]=backeuler(tspan,y0,h)
    if length(tspan)==2
        t0=tspan(1);
        tn=tspan(2);
    else
        error('MATLAB: Wrong Dimension Of Tspan');
    end
    n=floor((tn-t0)/h); 
    t=zeros(1,n+1);
    t(1)=t0; 
    y(:,1)=y0;
    for i=1:n
        t(i+1)=t(i)+h;
        y(:,i+1)=1/(1+h/t(i+1))*(y(:,i)+h*t(i+1)^3);
    end   
end

%% Heun method
function [t,y]=heun(ufunc,tspan,y0,h)
    if nargin<4
        h=0.01;
    end
    if length(tspan)==2
        t0=tspan(1);
        tn=tspan(2);
    else
        error('MATLAB: Wrong Dimension Of Tspan');
    end
    n=floor((tn-t0)/h); 
    t=zeros(1,n+1);
    t(1)=t0; 
    y(:,1)=y0;
    for i=1:n
        t(i+1)=t(i)+h;
        k1=ufunc(t(i),y(:,i));
        k2=ufunc(t(i)+h,y(:,i)+h*k1);
        y(:,i+1)=y(:,i)+h/2*(k1+k2);
    end
end

%% Ralston method
function [t,y]=ralston(ufunc,tspan,y0,h)
    if nargin<4
        h=0.01;
    end
    if length(tspan)==2
        t0=tspan(1);
        tn=tspan(2);
    else
        error('MATLAB: Wrong Dimension Of Tspan');
    end
    n=floor((tn-t0)/h); 
    t=zeros(1,n+1);
    t(1)=t0; 
    y(:,1)=y0;
    for i=1:n
        t(i+1)=t(i)+h;
        k1=ufunc(t(i),y(:,i));
        k2=ufunc(t(i)+h/2,y(:,i)+h/2*k1);
        k3=ufunc(t(i)+h*3/4,y(:,i)+h*3/4*k2);
        y(:,i+1)=y(:,i)+h/9*(2*k1+3*k2+4*k3);
    end
end

%% Runge-Kutta 4 method
function [t,y]=runge_kutta_4(ufunc,tspan,y0,h)
    if nargin<4
        h=0.01;
    end
    if length(tspan)==2
        t0=tspan(1);
        tn=tspan(2);
    else
        error('MATLAB: Wrong Dimension Of Tspan');
    end
    n=floor((tn-t0)/h); 
    t=zeros(1,n+1);
    t(1)=t0; 
    y(:,1)=y0;
    for i=1:n
        t(i+1)=t(i)+h;
        k1=ufunc(t(i),y(:,i));
        k2=ufunc(t(i)+h/2,y(:,i)+h/2*k1);
        k3=ufunc(t(i)+h/2,y(:,i)+h/2*k2);
        k4=ufunc(t(i)+h,y(:,i)+h*k3);
        y(:,i+1)=y(:,i)+h/6*(k1+2*k2+2*k3+k4);
    end
end