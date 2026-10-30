%{
%(1) Dinamica de modelo linear integrado numericamente
function [f,grad] = objfunc(t,u,ut,ufinal,modelo_x0,modelo_u0,modelo_par)
%A sua funcao objetivo no Spaider contem apenas a linha de comando
%apontando um valor para f. Para o meu caso tive que convocar a ode15s
%aqui dentro.
%Imagino que voce convoca a integracao numerica no arquivo function_ode.m.
%Se for isso, basta apagar as linhas abaixo convocando a ode15s.
%Tive que colocar valores vazios para os diversos parametros da modelo_eq.m
%Para uso com a Spaider vai precisar colocar os nomes das variaveis de volta.

Nut=size(ut,1);

if nargout > 1
    %[x,G] = dynmdl(ut,u,ut,ufinal,modelo_x0,modelo_u0,modelo_par);
    [x,G] = dynmdl(t,u,ut,ufinal,modelo_x0,modelo_u0,modelo_par);
    if size(x,1) < size(t,1) %if size(x,1) < Nut
        error('there was a failure in model integration');
    end
    %{
    %Gradient
    grad=zeros(Nut-1,1);
    for k=1:Nut-1
        idut=max(find(t<=ut(k)));
        idutp=max(find(t<=ut(k+1)));
        %grad(k)=sum(G(end,1,idut+1:idutp));
        grad(k)=sum(G(end,1,idut:idutp-1));
    end
    grad=grad*1e3 %grad=grad'*1e3; %Retornando ordem de grandeza
    %}
    %Gradient
    %grad=G(end,1,2:Nut);
    %grad=grad(:)'*1e3; %Retornando ordem de grandeza    
else
    x = modelo_eq(t,u,ut,ufinal,modelo_x0,modelo_u0,modelo_par);
    if size(x,1) < size(t,1) %if size(x,1) < Nut
        error('there was a failure in model integration');
    end
end
f=x(end,end);
end
%}
%{
%(2) Dinâmica de modelo linear via solução analítica linear impulso
function [f,grad] = objfunc(t,u,ut,ufinal,modelo_x0,modelo_u0,modelo_dFdx,modelo_dFdu)
%Funciona apenas quando tstep<=1 no arquivo main

Nut=size(ut,1);
[x,G] = dynmdl(u,u,ut,ufinal,modelo_x0,modelo_u0,modelo_dFdx,modelo_dFdu);
f = x(end,end);

%Gradient
grad=zeros(Nut-1,1);
for k=1:Nut-1
	idut=max(find(t<=ut(k)));
	idutp=max(find(t<=ut(k+1)));
	for j=idut:idutp-1
        grad(k) = grad(k) + G(end,1,j);
    end
end
grad=grad'*1e3; %Retornando ordem de grandeza
end
%}

%(3) Dinâmica de modelo linear via solução analítica linear degrau
function [f,grad] = objfunc(u,ut,ufinal,modelo_x0,modelo_u0,modelo_dFdx,modelo_dFdu)
Nut=size(ut,1);
[x,G] = dynmdl(ut,u,ut,ufinal,modelo_x0,modelo_u0,modelo_dFdx,modelo_dFdu);
%f = -x(end,end);
f = x(end,end);

%Gradient
%grad=-G(end,1,2:Nut);
grad=G(end,1,2:Nut);
grad=grad(:)'*1e3; %Retornando ordem de grandeza
end
%}