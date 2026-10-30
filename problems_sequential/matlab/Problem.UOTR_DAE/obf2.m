function [J dJ] = obf2(f,t_ult,t_out,df,ws,y,stages)
%global x
Nut=size(ut,1);
[x,G] = modelo_eq(t,u,ut,ufinal,modelo_x0,modelo_u0,modelo_dFdx,modelo_dFdu);
f = x(end,end);

%Gradient
dJ=zeros(Nut-1,1);
for k=1:Nut-1
	idut=max(find(t<=ut(k)));
	idutp=max(find(t<=ut(k+1)));
	for j=idut:idutp-1
        dJ(k) = dJ(k) + G(end,1,j);
    end
end
dJ=dJ'*1e3; %Retornando ordem de grandeza

%grad=G(end,1,1:Nut-1)*1e3; %Retornando ordem de grandeza
%grad=grad(:)';
end