%% Simulation of bouncing ball
% Author        Date        Version    Modification
% Liu Zhaoyong  2024-3-18   1.0        generate gif

%%
clc;clear;close all;

t0 = 0;
tf = 14;
y0 = [10; 0];
refine = 4;
options = odeset('Events',@events,'OutputSel',1,...
   'Refine',refine);

fig = figure("Color","white");
ax = axes;
ax.XLim = [t0 tf];
ax.YLim = [0 12];
box on
hold on;

tout = t0;
yout = y0';
teout = [];
yeout = [];
ieout = [];

collisionTimes = 12;
nImages = collisionTimes;

for i = 1:collisionTimes
   % Solve until the first terminal event.
   [t,y,te,ye,ie] = ode45(@f,[t0 tf],y0,options);
   
   if ~ishold
      hold on
   end
   % Accumulate output.  This could be passed out as output arguments.
   nt = length(t);
   tout = [tout; t(2:nt)];
   yout = [yout; y(2:nt,:)];
   teout = [teout; te];          % Events at tstart are never reported.
   yeout = [yeout; ye];
   ieout = [ieout; ie];
   
   plot(t,y(:,1),'b-o')
   if i ~= 1
    plot(t(1),y(1,1),'ro')
   end
   plot(te,ye(1),'ro')
   xlabel('Time (s)');
   ylabel('Height (m)');
   title('Ball trajectory and the events');
   hold off;

   frame = getframe(fig);
   im{i} = frame2im(frame);
   
    
   % Set the new initial conditions, with .9 attenuation.
   y0(1) = 0;
   y0(2) = -0.9*y(nt,2);
   
   % A good guess of a valid first timestep is the length of the last valid
   % timestep, so use it for faster computation.  'refine' is 4 by default.
   options = odeset(options,'InitialStep',t(nt)-t(nt-refine),...
      'MaxStep',t(nt)-t(1));
   
   t0 = t(nt);
end



filename = 'BouncingBall.gif'; % Specify the output file name
for j = 1:nImages
    [A,map] = rgb2ind(im{j},256);
    if j == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.3);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.3);
    end
end

%%
function dydt = f(t,y)
  dydt = [y(2); -9.8-0.5*0.05*sign(y(2))*y(2)^2];
end

%%
function [value,isterminal,direction] = events(t,y)
% Locate the time when height passes through zero in a decreasing direction
% and stop integration.
value = y(1);     % detect height = 0
isterminal = 1;   % stop the integration
direction = -1;   % negative direction
end