function dx = SLS_fun(t,x)
% 功能：dx = SLS_fun(t,x) 切换系统状态方程
% 输入：时间t,状态x
% 输出：导数dx
A1=[-1 10;-100 -1];
A2=[-1 100;-10 -1];
if x(1)*x(2)>=0
     A=A1;
else
     A=A2;
end
dx=A*x;