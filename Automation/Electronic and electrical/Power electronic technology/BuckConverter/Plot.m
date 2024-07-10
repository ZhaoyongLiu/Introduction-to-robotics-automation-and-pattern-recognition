% Plot the figures

t = out.tout;
uGE = out.uGE;
i = out.i;
v = out.v;
tinf = max(t);

figure(1)
subplot(2,1,1)
stairs(t,uGE,'Color','k','Linewidth',1.2)
xlabel('Time (s)')
ylabel('Voltage (V)')
title('uGE')
axis([-inf tinf -0.2 1.2])
subplot(2,1,2)
stairs(t(end-100:end),uGE(end-100:end),'Color','k','Linewidth',1.2)
xlabel('Time (s)')
ylabel('Voltage (V)')
title('uGE')
axis([-inf tinf -0.2 1.2])

figure(2)
subplot(2,1,1)
plot(t,i,'Color','magenta','Linewidth',1.2)
hold on;
ymax = 30.19;
ymin = 29.81;
yline([ymax ymin],'--',{'$I_{Max}$','$I_{Min}$'},'LineWidth',1.2,'Interpreter','Latex')
hold off;
xlabel('Time (s)')
ylabel('Current (A)')
title('Load Current')
axis([-inf tinf 29.5 30.5])
subplot(2,1,2)
plot(t(end-100:end),i(end-100:end),'Color','magenta','Linewidth',1.2)
hold on;
ymax = 30.19;
ymin = 29.81;
yline([ymax ymin],'--',{'$I_{Max}$','$I_{Min}$'},'LineWidth',1.2,'Interpreter','Latex')
hold off;
xlabel('Time (s)')
ylabel('Current (A)')
title('Load Current')
axis([-inf tinf 29.5 30.5])

figure(3)
subplot(2,1,1)
plot(t,v,'Color','b','Linewidth',1.2)
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Load Voltage')
axis([-inf tinf -10 110])
subplot(2,1,2)
plot(t(end-100:end),v(end-100:end),'Color','b','Linewidth',1.2)
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Load Voltage')
axis([-inf tinf -10 110])