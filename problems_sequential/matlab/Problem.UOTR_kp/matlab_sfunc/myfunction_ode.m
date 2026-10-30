%{
%Nao-linear
function [t_out,f_out,dFdx] = myfunction_ode(t_out,ws,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg,dFdx)
global modelo_x0 modelo_u0 modelo_par
%global modelo_dFdx
global JPatt

f_out=[]; %Initialize f_out

%States' initial conditions
%uf(1,:)=[ws(1)*1e3;modelo_u0(2:end)]';
x(1,:)=modelo_x0';
dx(1,:)=0*modelo_x0';

%Save to f_out
f_out(1,:)=x(end,:)';

%{
fhand=@(t,x)uotr_ctrl(t,x,uf(1,:)',1,[],modelo_par);
fval=fhand([],modelo_x0);
[dFdx_lin(:,:,1),~]=numjac(fhand,[],modelo_x0,fval,zeros(Nx,1));
modelo_dFdx=dFdx_lin(:,:,1);

fhand=@(t,u)uotr_ctrl(t,modelo_x0,[u;modelo_u0(2:end)],1,[],modelo_par);
fval=fhand([],modelo_u0(1));
[dFdu_lin(:,:,1),~]=numjac(fhand,[],modelo_u0(1),fval,zeros(Nu,1));
modelo_dFdu=dFdu_lin(:,:,1);
%}

%ode15s options and function handles
opt=odeset('RelTol',1e-5,'AbsTol',1e-6,'Stats','off','JPattern',JPatt);
%opt=odeset('RelTol',1e-5,'AbsTol',1e-6,'Stats','off','jacobian',Jac);
%opt=odeset('RelTol',1e-5,'AbsTol',1e-6,'Stats','off');
odefunc=@(t,x)modelo_eq(t,x,ws,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg);

%ODE integration
[~,f_out]=ode15s(odefunc,t_out,f_out(end,:)',opt);

%Create code to update jacobians
%{
for k=2:length(tspan)
    %Solve state + sensitivity ODE system from t(k-1) to t(k)
    opt=odeset('RelTol',1e-5,'AbsTol',1e-6,'Stats','off','Jacobian',Jac);
    odefunc=@(t,x)( [dxk(1:end-1);0;dsk] + modelo_eq(t,x,ws,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg) );
    [tv,xv]=ode15s(odefunc,[tspan(k-1) tspan(k)],[xlin(k-1,:)';s(k-1,:)'],opt);
    
    %Allocate t(k), xlin(k), s(k)
    t(k)=tv(end);
    xlin(k,1:Nx)=real(xv(end,1:Nx));
    s(k,:)=real(xv(end,Nx+1:end));
    
    %Update dFdx and dFdu
    dx_lin(k,:)=(xv(end,1:Nx)'-xv(end-1,1:Nx)')/(tv(end)-tv(end-1));
    dx_lin(k-1,:)=(xv(2,1:Nx)'-xv(1,1:Nx)')/(tv(2)-tv(1));
    uf(k,2:end)=modelo_u0(2:end)'; uf(k,1)=interpl(t(k),ws,ts,1e-2)*1e3;
    if abs((uf(k,:)-uf(k-1,:))./uf(k-1,:))>1e-3     %for abrupt change in u
        %Update dFdx
        fhand=@(t,x)uotr_ctrl(t,x,uf(k,:)',1,[],modelo_par);
        fval=fhand([],xlin(k,:)');
        [dFdx_lin(:,:,k),~]=numjac(fhand,[],xlin(k,:)',fval,zeros(Nx,1));
        if imag(dFdx_lin(:,:,k))>=1e-6
            error('Solution of states is imaginary')
        end
        dFdx_lin(:,:,k)=real(dFdx_lin(:,:,k)');
        modelo_dFdx=dFdx_lin(:,:,k);
        modelo_x0=xlin(k,:)';
        [T,D]=eig(modelo_dFdx); dzk=(T\dx_lin(k,:)');
        for j=1:Nx
            if real(D(j,j))>=-eps
                dzk(j)=0;
            end
        end
        dxk=T*dzk; dxk=real(dxk);
        
        %Update ODE system states + sensitivity jacobian
        Jac=modelo_dFdx;
        for j=1:n_c
            for is = 1:ns_a(j)
                Jac=blkdiag(Jac,modelo_dFdx);
            end
        end
        
        %Update dFdu
        fhand=@(t,u)uotr_ctrl(t,xlin(k,:)',[u;modelo_u0(2:end)],1,[],modelo_par);
        fval=fhand([],uf(k,1));
        [dFdu_lin(:,:,k),~]=numjac(fhand,[],uf(k,1),fval,zeros(Nu,1));
        dFdu_lin(:,:,k)=real(dFdu_lin(:,:,k));
        modelo_dFdu=dFdu_lin(:,:,k);
        modelo_u0=uf(k,:)';
    else
        if abs((dx_lin(k,:)'-dx_lin(k-1,:)')./dx_lin(k,:)')>1e-3    %for abrupt change in dx
            %Update dFdx
            fhand=@(t,x)uotr_ctrl(t,x,uf(k,:)',1,[],modelo_par);
            fval=fhand([],xlin(k,:)');
            [dFdx_lin(:,:,k),~]=numjac(fhand,[],xlin(k,:)',fval,zeros(Nx,1));
            if imag(dFdx_lin(:,:,k))>=1e-6
                error('Solution of states is imaginary')
            end            
            dFdx_lin(:,:,k)=real(dFdx_lin(:,:,k));
            modelo_dFdx=dFdx_lin(:,:,k);
            modelo_x0=xlin(k,:)';
            [T,D]=eig(modelo_dFdx); dzk=(T\dx_lin(k,:)');
            for j=1:Nx
                if real(D(j,j))>=-eps
                    dzk(j)=0;
                end
            end
            dxk=T*dzk; dxk=real(dxk);
            
            %Update ODE system states + sensitivity jacobian
            Jac=modelo_dFdx;
            for j=1:n_c
                for is = 1:ns_a(j)
                    Jac=blkdiag(Jac,modelo_dFdx);
                end
            end            
            
            %Update dFdu
            fhand=@(t,u)uotr_ctrl(t,xlin(k,:)',[u;modelo_u0(2:end)],1,[],modelo_par);
            fval=fhand([],uf(k,1));
            [dFdu_lin(:,:,k),~]=numjac(fhand,[],uf(k,1),fval,zeros(Nu,1));
            dFdu_lin(:,:,k)=real(dFdu_lin(:,:,k));
            modelo_dFdu=dFdu_lin(:,:,k);
            modelo_u0=uf(k,:)';
        else
            %No updates
            dFdx_lin(:,:,k)=dFdx_lin(:,:,k-1);
            dFdu_lin(:,:,k)=dFdu_lin(:,:,k-1);
        end
    end    
end
%}

%Sensitivity analytical solution
s=modelo_eqsens(t_out,ws,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg,dFdx);

%Stalling sensitivities into block vector in f_out
f_out=[f_out,s];
end
%}

%Linear
function [t_out,f_out] = myfunction_ode(t_out,ws,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg)
global modelo_x0 modelo_u0 modelo_dFdx JPatt

%Nx=size(modelo_x0,1); Nu=size(modelo_u0(1),1);

f_out=[]; %Initialize f_out

%States' initial conditions
%uf(1,:)=[ws(1)*1e3;modelo_u0(2:end)]';
x(1,:)=modelo_x0';

%{
%Sensitivity functions initial conditions
sv=[];
for j=1:n_c
    Id=zeros(n_s,ns_a(j));
    sv=[sv';Id(:)]';
end
s(1,:)=sv;
%}
%Save to f_out
f_out(1,:)=x(end,:)';
%f_out(1,:)=[x(end,:)';s(end,:)'];

%ODE system states + sensitivity jacobian
Jac=modelo_dFdx;
%{
for j=1:n_c
    for is = 1:ns_a(j)
        %Jac=blkdiag(Jac,dFdx);
        Jac=blkdiag(Jac,modelo_dFdx);
    end
end
%}
%ode15s options and function handles
%opt=odeset('RelTol',1e-5,'AbsTol',1e-6,'Stats','off','JPattern',JPatt);
opt=odeset('RelTol',1e-5,'AbsTol',1e-6,'Stats','off','jacobian',Jac);
%opt=odeset('RelTol',1e-5,'AbsTol',1e-6,'Stats','off');
odefunc=@(t,x)modelo_eq(t,x,ws,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg);

%ODE integration
[~,f_out]=ode15s(odefunc,t_out,f_out(end,:)',opt);

%Sensitivity analytical solution
s=modelo_eqsens(t_out,ws,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg);

%Stalling sensitivities into block vector in f_out
f_out=[f_out,s];
end
%}