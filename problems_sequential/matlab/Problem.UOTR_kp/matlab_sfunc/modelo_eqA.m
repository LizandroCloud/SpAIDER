%{
%Nao-linear
function [dx] = modelo_eq(t_out,x,ws,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg)
global modelo_u0 modelo_par
global modelo_dFdx modelo_dFdu

%Input interpolation
%[u] = ver_reg_est1(ws,yi,t_out);
%[uref] = ver_reg_est1(ws,yi,t_out);
uf=modelo_u0; uf(1)=interpl(t_out,ws,ts,1e-2)*1e3;

%States' derivative
dx=uotr_ctrl(t_out,x(1:n_s),uf,1,[],modelo_par);
dx=dx*tal;

%{
%{
%Compute jacobians for sensitivity matrices
fhand=@(t,x)uotr_ctrl(t,x,uf,1,[],modelo_par);
fval=fhand([],x(1:n_s));
if imag(fval)>=1e-6
    error('Solution of states is imaginary')
end                    
[dFdx,~]=numjac(fhand,[],x(1:n_s),real(fval),zeros(n_s,1));
if imag(dFdx)>=1e-6
    error('Solution of states is imaginary')
end            
dFdx=real(dFdx);

fhand=@(t,u)uotr_ctrl(t,x(1:n_s),[u;modelo_u0(2:end)],1,[],modelo_par);
fval=fhand([],uf(1));
if imag(fval)>=1e-6
    error('Solution of states is imaginary')
end                            
[dFdu,~]=numjac(fhand,[],uf(1),real(fval),zeros(n_c,1));
dFdu=real(dFdu);
%}
%Sensitivity matrices allocation in output [dx]
df=[modelo_dFdx modelo_dFdu*1e3]; %df=[dFdx dFdu*1e3];
jk=double(df);
j0=0;
jk(:,1:n_s)=df(:,1:n_s); %Jacobian dFdx
jk(:,n_s+1:n_s+n_c)=0;
for j=1:n_c
	for is = 1:ns_a(j)
        id=zeros(ns_a(j),1); id(is)=1;
		uws=interpl(t_out,id,ts,1e-2);
		jk(:,n_s+j)=uws*df(:,n_s+j); %jk(:,n_s+j)=uws(j,is)*df(:,n_s+j);
		k=1:(n_s);
		dx(n_s*is+k+j0) = jk(k,1:n_s)*x(n_s*is+1+j0:n_s*is+n_s+j0) + jk(k,n_s+j);	%Sensitivity dsdt=dFdx'*s+dFdu';
        dx(n_s*is+k+j0)=dx(n_s*is+k+j0)*tal;
	end 
	j0=length(dx)-n_s;
end
%}
end
%}

%Linear
function [dx] = modelo_eq(t_out,x,ws,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg)
global modelo_x0 modelo_u0 modelo_dFdx modelo_dFdu

%Nx=size(modelo_x0,1);	%Number of states

%Input interpolation
%[u] = ver_reg_est1(ws,yi,t_out);
%[uref] = ver_reg_est1(ws,yi,t_out);
[u] = interpl(t_out,ws,ts,1e-2);

%Model and jacobian
%{
dx=modelo_dFdx*(x-modelo_x0)+modelo_dFdy*(u*1e-3-modelo_u0(1));
dx=(T*D/T)*(x-modelo_x0)+modelo_dFdy*(u*1e-3-modelo_u0(1));
(T\dx)=eig(modelo_dFdx)*(T\(x-modelo_x0))+T\modelo_dFdy*(u*1e-3-modelo_u0(1));
dz=D*(z-z0)+B*(u*1e-3-modelo_u0(1));
%}

[T,D]=eig(modelo_dFdx); B=T\modelo_dFdu;
z=T\x(1:n_s); z0=T\modelo_x0(1:n_s);
for j=1:n_s
	if real(D(j,j))>=-1e-6
		dz(j)=B(j,:)*u*1e3;
		dfz(j,:)=[zeros(1,n_s) B(j,:)*1e3];
	else
		dz(j)=D(j,j)*(z(j)-z0(j))+B(j,:)*(u*1e3-modelo_u0(1));
		Id=zeros(1,n_s); Id(j)=D(j,j);
		dfz(j,:)=[Id B(j,:)*1e3];
	end
end
dz=dz(:); dx=T*dz*tal;
if imag(dx)>=1e-6
    error('Solution of states is imaginary')
end
dx=real(dx);
%{
%Sensitivity matrices allocation in output [dx]
df=[T*dfz(:,1:Nx)/T T*dfz(:,Nx+1:end)]; df=real(df);
jk=double(df);
j0=0;
jk(:,1:n_s)=df(:,1:n_s); %Jacobian dfdx
jk(:,n_s+1:n_s+n_c)=0;
for j=1:n_c
	for is = 1:ns_a(j)
        id=zeros(ns_a(j),1); id(is)=1;
		uws=interpl(t_out,id,ts,1e-2);
		jk(:,n_s+j)=uws*df(:,n_s+j); %jk(:,n_s+j)=uws(j,is)*df(:,n_s+j);
		k=1:(n_s);
		dx(n_s*is+k+j0) = jk(k,1:n_s)*x(n_s*is+1+j0:n_s*is+n_s+j0) + jk(k,n_s+j);	%Sensitivity dsdt=dfdx'*s+dfdu';
        dx(n_s*is+k+j0)=dx(n_s*is+k+j0)*tal;
	end 
	j0=length(dx)-n_s;
end
%}
end
%}