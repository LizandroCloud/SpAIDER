function [ dx ] = modelo_eq_sens( t_out, x, ws, ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg)
global modelo_x0 modelo_u0 modelo_dFdx modelo_dFdu

Nx=size(modelo_x0,1);	%Number of states

%Input interpolation
[u] = ver_reg_est1(ws,yi,t_out);
%[uref] = ver_reg_est1(ws,yi,t_out);
% [u] = interp1(ts,ws,t_out);


%Model and jacobian
%{
dx=modelo_dFdx*(x-modelo_x0)+modelo_dFdy*(u*1e-3-modelo_u0(1));
dx=(T*D/T)*(x-modelo_x0)+modelo_dFdy*(u*1e-3-modelo_u0(1));
(T\dx)=eig(modelo_dFdx)*(T\(x-modelo_x0))+T\modelo_dFdy*(u*1e-3-modelo_u0(1));
dz=D*(z-z0)+B*(u*1e-3-modelo_u0(1));
%}

[T,D]=eig(modelo_dFdx); B=T\modelo_dFdu;
z=T\x(1:Nx); z0=T\modelo_x0(1:Nx);
for j=1:Nx
	if D(j,j)>=-1e-6
		dz(j)=B(j,:)*u*1e3;
		dfz(j,:)=[zeros(1,Nx) B(j,:)*1e3];
	else
		dz(j)=D(j,j)*(z(j)-z0(j))+B(j,:)*(u*1e3-modelo_u0(1));
		Id=zeros(1,Nx); Id(j)=D(j,j);
		dfz(j,:)=[Id B(j,:)*1e3];
	end
end
dz=dz(:);
dx=T*dz*tal;
df=[T*dfz(:,1:Nx)/T T*dfz(:,Nx+1:end)]*tal;
if imag(dx)>=1e-6
    error('Solution ofr states is imaginary')
else
   dx=real(dx);
   df=real(df);
end

%Sensitivity matrices allocation in output [dx]
jk=double(df);
j0=0;
jk(:,1:n_s)=df(:,1:n_s); %Jacobian dfdx
jk(:,n_s+1:n_s+n_c)=0;
for j=1:n_c
	for is = 1:ns_a(j)+1
        id=zeros(ns_a(j),1); id(is)=1;
		% uws=ver_reg_est1(ws,yi,t_out);
		jk(:,n_s+j)=uws*df(:,n_s+j); %jk(:,n_s+j)=uws(j,is)*df(:,n_s+j);
		k=1:(n_s);
		dx(n_s*is+k+j0) = jk(k,1:n_s)*x(n_s*is+1+j0:n_s*is+n_s+j0) + jk(k,n_s+j);	%Sensitivity dsdt=dfdx'*s+dfdu';
	end 
	j0=length(dx)-n_s;
end
end