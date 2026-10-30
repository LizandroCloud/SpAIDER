function [dx] = modelo_eqX(t_out,x,ws,ks,yi)
global modelo_x0 modelo_u0 modelo_dFdx modelo_dFdu
                 
Nx=size(modelo_x0,1);	%Number of states
tal=1;
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
if imag(dx)>=1e-6
    error('Solution ofr states is imaginary')
else
   dx=real(dx);
end
end