clc; clear; close all;
x = 0:2*pi;

y = sin(x);
xx = 0.5:2*pi
yTrue = sin(xx)
yLinear = interp1(x,y,xx)
yLag = Laginterp(x,y,xx)
L = Laginterp(x,y);


figure(1)
plot(x,y,'s','MarkerSize',10,'MarkerEdgeColor','b','MarkerFaceColor',[0.5,0.5,0.5]) 
hold on
fplot(@(x) sin(x), [0,2*pi],'k','LineWidth',1.2)
fplot(L, [0,2*pi],'r','LineWidth',1.2)
legend('Interpolation nodes','$y = \sin x$','$y = L(x)$','Interpreter','latex','Fontsize',15)
set(gca,'Fontsize', 12)
axis([0 2*pi -inf inf])


function L = Laginterp(x1, y1, x0)
% Function: Calculate Lagrange polynomial interpolation function
% Input:  x1  -- interpolation nodes
%         y1  -- function values y1 = f(x1)
%         x0  -- nodes to be interpolated
% Output: L  -- Lagrange polynomial interpolation function L(x) = y0*l0(x)+ ... + yk*lk(x)
% Syntax: L = Laginterp(x1, y1) -- Given interpolation nodes x and function values y = f(x),
%                                calculate the Lagrange interpolation function L(x)
%         L = Laginterp(x1, y1, x0) -- Calculate L(x) and the interpolation values at x0
nx = length(x1);
ny = length(y1);
if nx ~= ny
   warning('The length of x1 and y1 must be the same')
   return;
end
if nargin == 2
   syms x;
   w = 0;
   digits(6)
   for k = 1:nx
       u = 1.0;
       for j = 1:nx
           if j ~= k
               u = u*(x-x1(j))/(x1(k)-x1(j));
           end
       end
       w = w + u*y1(k);
   end
   L = simplify(w);
   disp('拉格朗日插值函数的小数形式为')
   vpa(L)
else
    m = length(x0);
    for i = 1:m
        t = 0.0;
        for k = 1:nx
            u = 1.0;
            for j = 1:nx
                if j~= k
                    u = u*(x0(i)-x1(j))/(x1(k)-x1(j));
                end
            end
            t = t + u*y1(k);
        end
        L(i) = t;
    end

end  

end
