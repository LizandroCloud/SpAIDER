%{
%(1) Dinâmica de modelo linear diferencial para integração numérica
function dx = modelo_eq(t_out,x_out,ws,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,tal,est,yi,utarg,tstarg,nstarg)
global modelo_dFdx modelo_dFdu modelo_x0 modelo_u0 ufinal

%Comentei a sua regfunc e to usando a minha porque não tenho o arqui-
%vo .m da sua.
%OBS: Perceba que estou usando [ws,ufinal] como argumento da regfunc,
%onde ufinal é um valor terminal fixo para a vazão de entrada. 
ws_reg=interpl(t_out,[ws;ufinal],1e-2); %ws_reg=ver_reg_est1([ws;ufinal],yi,t_out);

u_out=modelo_u0;
u_out(1)=ws_reg*1e3; %Retornando ordem de grandeza de u para uso correto no modelo.
dx=modelo_dFdx*(x_out-modelo_x0)+modelo_dFdu*(u_out-modelo_u0);
dx=dx*tal;

end
%}

%(2) Dinâmica de modelo linear via solução analítica
%function [x,grad] = modelo_eq(u)
function [x,G] = modelo_eq(t,u,ut,ufinal,modelo_x0,modelo_u0,modelo_dFdx,modelo_dFdu)
%global t ut x ufinal modelo_x0 modelo_u0 modelo_dFdx modelo_dFdu
Nx=size(modelo_x0,1); Nu=size(modelo_u0,1); kspan=size(t,1); %Nut=size(ut,1);

%Eigenvalues and eigenvectors
%[T,D]=eig(modelo_dFdx); d=diag(D);
%id=find(real(d)<=-eps); %idi=find(real(d)>-eps); Nxd=size(id,1);

%Control horizon dynamics and transition matrix
U=ones(kspan,1)*modelo_u0';
for k=1:kspan
	Phi(:,:,k)=expm(modelo_dFdx*(t(k)-t(1)));
	G(:,:,k)=Phi(:,:,k)*modelo_dFdu;
    %Phi(:,:,k)=expm(modelo_dFdx(1:Nx-1,1:Nx-1)*(t(k)-t(1)));
	%Phi(:,:,k)=expm(D(id,id)*(t(k)-t(1)));
	%U(k,1)=interpl(t(k),[u;ufinal],ut,1e-2)*1e3; %for guess 1. Retornando ordem de grandeza de u para uso correto no modelo.
	U(k,1)=interpl(t(k),[u;u(end)],ut,1e-2)*1e3; %for guess 2. Retornando ordem de grandeza de u para uso correto no modelo.
end

%Convolution
cnvlt=zeros(kspan,Nx);
%zcnvlt=zeros(Nxt,Nx-1);
for k=2:kspan
	%if k==1
	%	zcnvlt(k,:) = (-D(id,id)\(eye(Nx-1)-Phi(:,:,k))*(T(id,id)\modelo_dFdu(1:Nx-1,:))*(U(k,:)'-modelo_u0))';
	%else
		%cnvlt(k,:) = (-(T*inv(D)/T)*((T*Phi(:,:,k-1)/T)-(T*Phi(:,:,k)/T))*modelo_dFdu*(U(k,:)'-modelo_u0))';
		%T\cnvlt(k,:) = (-inv(D)/T)*((T*Phi(:,:,k-1)/T)-(T*Phi(:,:,k)/T))*modelo_dFdu*(U(k,:)'-modelo_u0))';
	%	zcnvlt(k,:) = (-D(id,id)\(Phi(:,:,k-1)-Phi(:,:,k))*(T(id,id)\modelo_dFdu(1:Nx-1,:))*(U(k,:)'-modelo_u0))';
        %zcnvlt(k,:) = (-D(id,id)\(eye(Nx-1)-Phi(:,:,k))*(T(id,id)\modelo_dFdu(1:Nx-1,:))*(U(k,:)'-modelo_u0))';
	%end
    %zcnvlt(k,:) = (-D(id,id)\(Phi(:,:,k-1)-Phi(:,:,k))*(T(id,id)\modelo_dFdu(1:Nx-1,:))*(U(k,:)'-modelo_u0))';
    %zcnvlt(k,:) = (-D(id,id)\(eye(Nx-1)-Phi(:,:,k))*(T(id,id)\modelo_dFdu(1:Nx-1,:))*(U(k,:)'-modelo_u0))';
    %zcnvlt(k,:) = zcnvlt(k-1,:) + (Phi(:,:,k)*(T(id,id)\modelo_dFdu(1:Nx-1,:))*(U(k,:)'-modelo_u0))';
    %cnvlt(k,:) = cnvlt(k-1,:) + (Phi(:,:,k)*modelo_dFdu(1:Nx-1,:)*(U(k,:)'-modelo_u0))';
	%cnvlt(k,:) = cnvlt(k-1,:) + (Phi(:,:,k)*modelo_dFdu*(U(k,:)'-modelo_u0))';
	cnvlt(k,:) = cnvlt(k-1,:) + (G(:,:,k)*(U(k,:)'-modelo_u0))';
end

%States
%z=ones(Nxt,1)*((T(id,id)\modelo_x0(1:Nx-1)))';
x=ones(kspan,1)*modelo_x0';
for k=1:kspan
	%x(k+1,:) = (Phi(:,:,k+1)*(x(1,:)'-modelo_x0) + cnvlt(k,:)')';
    %z(k,:) = z(1,:) + (Phi(:,:,k)*(z(1,:)'-T(id,id)*modelo_x0(1:Nx-1)) + zcnvlt(k,:)')';
	%x(k,1:Nx-1) = modelo_x0(1:Nx-1)' + real((T(id,id)\(Phi(:,:,k)*(z(1,:)'-T(id,id)*modelo_x0(1:Nx-1)) + zcnvlt(k,:)'))');
    %x(k,1:Nx-1) = modelo_x0(1:Nx-1)' + real((T(id,id)*zcnvlt(k,:)'))';
    %x(k,1:Nx-1) = modelo_x0(1:Nx-1)' + cnvlt(k,:);
	%x(k,end) = ((t(k)-t(1))*modelo_dFdu(end,:)*(U(k,:)'-modelo_u0))';
    x(k,:) = modelo_x0' + cnvlt(k,:);
end
%x(:,end) = -x(:,end);

%Gradient
%grad=zeros(Nut-1,1);
%G=zeros(Nx,Nu,Nut-1);
%for k=1:Nut-1
%	idut=max(find(t<=ut(k))); idutp=max(find(t<=ut(k+1)));
%	for j=idut:idutp-1
%        G(:,:,k) = G(:,:,k) + Phi(:,:,j)*modelo_dFdu;
%    end
	%grad(k)=G(end,1,k)*1e3; %Retornando ordem de grandeza
%end
%end

%{
%(3) Dinâmica de modelo não-linear via solução analítica linearizado por partes
function [x,grad] = modelo_eq(u)
global t ut x ufinal modelo_x0 modelo_u0 modelo_par
Nx=size(modelo_x0,1); Nu=size(modelo_u0,1); Nxt=size(t,1); Nut=size(ut,1);

U=ones(Nxt,1)*modelo_u0';
cnvlt=zeros(Nxt,Nx);
x(1,1:Nx)=modelo_x0';
Phi(1:Nx,1:Nx,1)=eye(Nx);
G(1:Nx,1:Nu,1)=zeros(Nx,Nu);
for k=1:Nxt-1
	%Input interpolation
	%U(k,1)=interpl(t(k),[u;ufinal],ut,1e-2)*1e3; %for guess 1. Retornando ordem de grandeza de u para uso correto no modelo.
	U(k,1)=interpl(t(k),[u;u(end)],ut,1e-2)*1e3; %for guess 2. Retornando ordem de grandeza de u para uso correto no modelo.

	%Jacobian dFdx 
	fhand=@(ts,xs)uotr_ctrl(ts,xs,U(k,:)',1,[],modelo_par);
	fval=fhand(0,x(k,:)');
	[modelo_dFdx(:,:,k),~]=numjac(fhand,0,x(k,:)',fval,zeros(Nx,1));

	%Jacobian dFdu
	fhand=@(ts,us)uotr_ctrl(ts,x(k,:)',us,1,[],modelo_par);
	fval=fhand(0,U(k,:)');
	[modelo_dFdu(:,:,k),~]=numjac(fhand,0,U(k,:)',fval,zeros(Nu,1));	

	%Transition matrix
	Phi(:,:,k+1)=expm(modelo_dFdx(:,:,k)*(t(k+1)-t(k)))*Phi(:,:,k);

	%Convolution
	for j=1:k
		G(:,:,k+1) = (t(j+1)-t(j))*Phi(:,:,k-j+1)*modelo_dFdu(:,:,j);
		cnvlt(k+1,:) = (cnvlt(k+1,:)' + G(:,:,k+1)*(U(j+1,:)'-U(j,:)'))';
	end

	%k=5:
	%cnvlt(k,:) = ...
	%	(t(2)-t(1))*Phi(:,:,4)*modelo_dFdu(:,:,1)*(U(2,:)'-U(1,:)')' + ...
	%	(t(3)-t(2))*Phi(:,:,3)*modelo_dFdu(:,:,2)*(U(3,:)'-U(2,:)')' + ...
	%	(t(4)-t(3))*Phi(:,:,2)*modelo_dFdu(:,:,3)*(U(4,:)'-U(3,:)')' + ...
	%	(t(5)-t(4))*Phi(:,:,1)*modelo_dFdu(:,:,4)*(U(5,:)'-U(4,:)')' + ...

	%States
	%x(k+1,:) = (Phi(:,:,k+1)*x(1,:)' + cnvlt(k+1,:)')';
	x(k+1,:) = modelo_x0' + cnvlt(k+1,:);
end

%Gradient
grad=zeros(Nut-1,1);
for k=1:Nut-1
	idt=max(find(t<=ut(k))); idtp=max(find(t<=ut(k+1)));
	dDUdu=zeros(Nxt-1,1); dDUdu(idt:idtp)=1;
	g=G(end,1,1:Nxt-1); g=g(:);
	grad(k)=g'*dDUdu*1e3; %Retornando ordem de grandeza
end
end
%}
%{
%(X) Dinâmica de modelo não-linear em malha aberta
function dx = modelo_eq(t_out,x,ws,ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,tal,est,yi,utarg,tstarg,nstarg)
global modelo_u0 modelo_par

ws_reg=ver_reg_est1(ws,yi,t_out);

%idx=[2,3,5,6,7,9,10,13,15,17,19,20,22,24,25,28];
%idy=[1,4,8,12,14,16,18,26]; %idy=[1,4,8,11,12,14,16,18,21,23,26,27];

idu=[1,5,7,9,15,16,19,20,25]; %idu=[1,5,7,9,14,15,16,19,20,21,23,25,26];
%idw=[2,3,4,6,8,10,11,12,13,17,18,22,24,27];

u=modelo_u0; u(idu)=ws_reg;
dx=uotr(t_out,x,u,1,[],modelo_par);
dx=dx*tal;

end
%}