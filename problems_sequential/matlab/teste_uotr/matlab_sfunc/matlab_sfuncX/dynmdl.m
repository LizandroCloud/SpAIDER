function x = dynmdl(u,uf)
global t Nd u0 x0 dFdx par
Nx=size(x0,1);
tstart=0; tpred=1*t(end); tstep=60; tspan=(tstart:tstep:tpred)';
opt=odeset('RelTol',1e-5,'AbsTol',1e-8,'Jacobian',dFdx);
[ts,x]=ode15s(@(ts,x)modelo_eq(ts,x,u,uf),tspan,x0,opt);
if ts(end) < tpred %|| tsb(end) < tpred
	error('there was a failure in model integration');
end
end

function dxdt = dynf(ts,xdyn,u,uf)
global dFdx dFdu x0 u0 par
Nu=size(u0,1);
ureg=interpl(ts,[u;uf],1e-2);
u2f=u0; u2f(1)=ureg*1e3;
dxdt=dFdx*(xdyn-x0)+dFdu*(u2f-u0);
end
