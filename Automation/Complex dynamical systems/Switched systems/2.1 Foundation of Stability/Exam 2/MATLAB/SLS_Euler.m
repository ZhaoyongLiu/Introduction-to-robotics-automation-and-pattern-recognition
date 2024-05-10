function [t,x,swi_t] = SLS_Euler(fun, tspan, x0, h)
% 功能：[t x] = SLS_Euler(fun,tspan,x0,h) 前向欧拉法求解常微分方程
% 输入：函数名fun,时间区间tspan,初值x0,步长h
% 输出：时间序列t,状态轨迹x,切换信号st
a=tspan(1);b=tspan(end);
M = floor(b-a)/h ;      % 离散点个数 M+1
x =zeros(length(x0),M+1);   % 行向量
t = a:h:b;
x(:,1) = x0;
swi_t=zeros(2,length(t)); % 切换信号
swi_t(1,:)=t; swi_t(2,1)=1;
for i = 1:M
    [dx,st]=fun(t(i), x(:,i));
    x(:,i+1) = x(:,i) +h *dx;
    swi_t(2,i+1)=st;
end