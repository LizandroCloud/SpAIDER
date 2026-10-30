function dydt=ufcvd(t,y)
global cA0 cB0 cAf;
%y(1) is V, y(2) is Vw, y(3) is cA, y(4) is cB
A=1.2; Vi=0.05; cBf=0.2; cAf=30;
k=5e-6;clim=280; sA=0.9;sB=0.15;V0=0.1;
J=k*log(clim/sA/y(3)-(1-sA)/sA);
Vw=Vi*(sA*log(cB0/cBf)-sB*log(cA0/cAf))/(sA-sB);
if and(y(1)<=Vi,y(2)<Vw);
alpha=1;
else
alpha=0;
end
Veqn=(alpha-1)*J*A;
Vweqn=alpha*J*A;
cAeqn=y(3)/y(1)*J*A*(sA-alpha);
cBeqn=y(4)/y(1)*J*A*(sB-alpha);
dydt=[Veqn;Vweqn;cAeqn;cBeqn];