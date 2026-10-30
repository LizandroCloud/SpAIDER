function [f,J,gradJ,c,ceq,gradc,gradceq] = computeall(t_out,ws,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg)
global modelo_x0

Nx=size(modelo_x0,1); 	%Number of states
Nt	= size(t_out,1);	%Number of time steps

%Solve dynamic model
funh_model=@(t,x)modelo_eq(t,x,ws,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg);
s0=[];
for j=1:n_c
    Id=zeros(n_s,ns_a(j));
    s0=[s0;Id(:)];
end
[t_out,f_out]=ode15s(funh_model,t_out,[modelo_x0;s0],odeset('RelTol',1e-5,'AbsTol',1e-8,'Stats','off'));

%Sensitivity matrix
j0=0;
for j=1:n_c
    for is = 1:ns_a(j)
		%k=1:(n_s);
		df(:,:,is+j0)=f_out(:,n_s*is+1:n_s*is+n_s);
    end
    j0=ns_a(j);
end

%States
f=f_out(:,1:Nx);

%Compute objective function and gradients
[J,gradJ] = fobj(f,[],t_out,df,[],[],[]);

%Compute constraints and gradients
[c,ceq,gradc,gradceq] = cr2(f,[],t_out,ts,df,ws,[],[]);
	
pause(1) % Simulate expensive computation
end