%{
%(1) Dinâmica de modelo linear integrado numericamente
function f = fobj(u)
global t x modelo_dFdx modelo_x0

%x=dynmdl(u,uf);

%A sua função objetivo no Spaider contém apenas a linha de comando
%apontando um valor para f. Para o meu caso tive que convocar a ode15s
%aqui dentro.
%Imagino que você convoca a integração numérica no arquivo function_ode.m.
%Se for isso, basta apagar as linhas abaixo convocando a ode15s.
%Tive que colocar valores varios para os diversos parâmetros da modelo_eq.m
%Para uso com a Spaider vai precisar colocar os nomes das variáveis de volta.

opt=odeset('RelTol',1e-5,'AbsTol',1e-8,'Jacobian',modelo_dFdx);
fhand=@(t_out,x_out)modelo_eq(t_out,x_out,u,[],[],[],[],[],[],[],[],[],[],[],1,[],[],[],[]);
[t_out,x]=ode15s(fhand,t,modelo_x0,opt);
if t_out(end) < t(end)
	error('there was a failure in model integration');
end

%OBS: Para o modelo linear, os valores de f serão negativos e sua magnitude precisa ser maximizada (min -f)
f=x(end,end);
end
%}
%{
%(2) Dinâmica de modelo linear via solução analítica
function [f,grad] = fobj(u)
global ut x ufinal modelo_x0 modelo_u0 modelo_dFdx modelo_dFdu
Nx=size(modelo_x0,1); Nu=size(modelo_u0,1); Nut=size(ut,1);

%Control horizon dynamics
U=ones(Nut,1)*modelo_u0';
U(1:Nut,1)=[u;ufinal]*1e3; %for guess 1. Retornando ordem de grandeza de u para uso correto no modelo.
%U(1:Nut,1)=[u;u(end)]*1e3; %for guess 2. Retornando ordem de grandeza de u para uso correto no modelo.

%Transition matrix
for k=1:Nut
	Phi(:,:,k)=expm(modelo_dFdx*(ut(k)-ut(1)));
end

%Convolution
cnvlt=zeros(Nx,1);
for k=1:Nut-1
	cnvlt = cnvlt + (ut(k+1)-ut(k))*Phi(:,:,Nut-k+1)*modelo_dFdu*(U(k,:)'-modelo_u0);
end

%Objective function
%x = Phi(:,:,end)*modelo_x0 + cnvlt;
x = modelo_x0 + cnvlt;
f = x(end);


%Gradient
grad=zeros(Nut-1,1);
for k=1:Nut-1
	G(:,:,k)=(ut(k+1)-ut(k))*Phi(:,:,Nut-k+1)*modelo_dFdu;
	grad(k)=G(end,1,k)*1e3; %Retornando ordem de grandeza
end
end
%}

%(3) Dinâmica de modelo não-linear via solução analítica linearizado por partes
function [f,grad] = fobj(u)
global x

[x,grad] = modelo_eq(u);
f = x(end,end);

end