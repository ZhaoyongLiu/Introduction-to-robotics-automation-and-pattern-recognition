%% Plot
t=out.tout;
y=out.simout.signals.values;

figure(1)
plot(t,y(:,1),'Color','k','Linewidth',1.2);
hold on;
plot(t,y(:,2),'Color','#77AC30','Linestyle','--','Linewidth',1.2);
legend('$x_1$','$x_2$','Interpreter','Latex','Fontsize',14);
xlabel('Time (s)');

figure(2)
f=arrowPlot(y(:,1),y(:,2),'number',1);
axis([-2 12 -4 10]);
xlabel('$x_1$','Interpreter','Latex','Fontsize',14);
ylabel('$x_2$','Interpreter','Latex','Fontsize',14);
set(f,'LineWidth',1.2)