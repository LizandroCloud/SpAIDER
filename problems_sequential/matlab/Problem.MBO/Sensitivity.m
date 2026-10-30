global wsx tex ksx nex
% k=1;
% wsx = [0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4];
% for t=0:1:10800
% [ux] = regx(wsx,t,tex,ksx,nex);
% ut(k)=ux;
% tx(k)=t;
% k=k+1;
% end
% plot(tx,ut);
nex=8;
tex = linspace(0,10800,8)
axc=[];
for i=1:10
    global wsx
    wsx(2:2) = 0.5*(1+0.015*(i-1));
    y0(210) = wsx(1);
    [tspan,zs] = ode15i('daesys_ss',[0:1:10800],y0,ydot0);
    axc(i)=zs(end,148)/0.001;
    b=length(zs(:,148));
    Jx=(zs(:,148)-ones(b,1)*0.005).^2 ;
    J(i)=trapz(Jx);
end

stairs(J)