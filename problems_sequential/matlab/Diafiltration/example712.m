function example712
global cA0 cB0
cA0=10; cB0=1.0;
options=odeset('events', @eventfcn, 'Maxstep',1);
[t,y]=ode15s(@ufcvd,[0 36000],[0.1 0 10 1], options);
cA=y(:,3)./cA0; cB=y(:,4)./cB0;
time=t/60;
timeout=t(end)/60
plot(time,cA,time,cB);
xlabel('Time (mins)');ylabel('Relative Conc. (g/L)')



