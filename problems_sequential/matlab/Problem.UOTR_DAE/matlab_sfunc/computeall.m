%{
%Nao-linear

function [f,J,gradJ,c,ceq,gradc,gradceq,dFdx] = computeall(t_out,ws,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg,GradObj,GradConstr,dFdx)
J=[]; gradJ=[]; c=[]; ceq=[]; gradc=[]; gradceq=[];

tStart=tic;
disp('Starting model integration...');
%ODE system: Compute states, states' derivatives, sensitivity matrices and update jacobians
[~,f_out,dFdx]=myfunction_ode(t_out,ws,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg,dFdx);
elaptime=toc(tStart);
disp(strjoin({'Integration finished succssesfully after ',num2str(elaptime,2),' seconds'},''));

%Get sensitivity matrices
j0=0;
for j=1:n_c
    for is = 1:ns_a(j)
		df(:,:,is+j0)=f_out(:,n_s*is+1:n_s*is+n_s);
    end
    j0=ns_a(j);
end

%Get states
f=f_out(:,1:n_s);

%{
%Compute ODE jacobian
if flagstart
    dx_old=dx;
    flagstart=false;
end
tic
[dFdx,dFdu]=myfunction_odejac(t_out,f,ws,ts,n_s,n_c,dx,dx_old,dFdx,dFdu);
elaptime=toc;
disp(strjoin({'Elapsed time for jacobian update = ',num2str(elaptime)},''));
%}

%Compute objective function and gradients
if GradObj 
    [J,gradJ] = obf2(f,[],t_out,df,[],[],[]);
else
    J = obf1(f,[],t_out,df,[],[],[]);
end

%Compute constraints and gradients
if GradConstr
    [c,ceq,gradc,gradceq] = cr2(f,[],t_out,ts,df,ws,[],[]);
else
    [c,ceq] = cr1(f,[],t_out,ts,df,ws,[],[]);
end
	
pause(1) % Simulate expensive computation
end
%}

%Linear
function [f,J,gradJ,c,ceq,gradc,gradceq] = computeall(t_out,ws,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg,GradObj,GradConstr)
J=[]; gradJ=[]; c=[]; ceq=[]; gradc=[]; gradceq=[];

tStart=tic;
disp('Starting model integration...');
%ODE system: Compute states, states' derivatives, sensitivity matrices and update jacobians
[~,f_out]=myfunction_ode(t_out,ws,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg);
elaptime=toc(tStart);
disp(strjoin({'Integration finished succssesfully after ',num2str(elaptime,2),' seconds'},''));

%Get sensitivity matrices
j0=0;
for j=1:n_c
    for is = 1:ns_a(j)
		df(:,:,is+j0)=f_out(:,n_s*is+1:n_s*is+n_s);
    end
    j0=ns_a(j);
end

%G.et states
f=f_out(:,1:n_s);

%Compute objective function and gradients
if GradObj 
    [J,gradJ] = obf2(f,[],t_out,df,[],[],[]);
else
    J = obf1(f,[],t_out,df,[],[],[]);
end

%Compute constraints and gradients
if GradConstr
    [c,ceq,gradc,gradceq] = cr2(f,[],t_out,ts,df,ws,[],[]);
else
    [c,ceq] = cr1(f,[],t_out,ts,df,ws,[],[]);
end
	
pause(1) % Simulate expensive computation
end
%}