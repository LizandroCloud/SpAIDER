function [c, ceq] = mycon_reactor61_haar(w, m, Hi, iDWT, p)
% Restricoes - Exemplo 6.1, versao Haar (Hi/iDWT)
%
% w = [tf; xdA(1:m); xdB(1:m); u(1:m)]
%
% Tempo normalizado tau em [0,1], t = tau*tf. Reconstrucao via Hi/iDWT
% (mesmo papel de myfun_vo.m / mycon_case63_haar.m):
%   cA(tau) = cAo + (xdA*inv(Hi))*iDWT
%   cB(tau) = cBo + (xdB*inv(Hi))*iDWT
%   V(tau)  = Vo  + tf*((u*inv(Hi))*iDWT)      <- V_dot=u e' linear,
%                                                  entao V vem direto da
%                                                  integral de u (nao
%                                                  precisa de uma
%                                                  variavel xdV separada)
%
% ceq: residuo da EDO para cA e cB (com o fator tf, pois d/dtau = tf*d/dt)
% c:   restricoes de caminho/terminal (desigualdade, <=0):
%        cB(tau) <= cBmax,  V(tau) <= Vmax,  nCdes <= nC(tf)

    tf  = w(1);
    xdA = w(2 : 1+m);
    xdB = w(2+m : 1+2*m);
    u   = w(2+2*m : 1+3*m);

    Hi_inv = inv(Hi); %#ok<MINV>

    cA = p.cAo + (xdA * Hi_inv) * iDWT;
    cB = p.cBo + (xdB * Hi_inv) * iDWT;
    V  = p.Vo  + tf * ((u * Hi_inv) * iDWT);

    dcA = -p.k*cA.*cB - (u./V).*cA;
    dcB = -p.k*cA.*cB + (u./V).*(p.cBin - cB);

    ceq = [xdA - tf*dcA, xdB - tf*dcB];

    nC_final = p.cAo*p.Vo - cA(end)*V(end);
    g_nC = p.nCdes - nC_final;   % <=0
    g_cB = cB - p.cBmax;         % <=0
    g_V  = V  - p.Vmax;          % <=0

    c = [g_nC, g_cB, g_V];

end
