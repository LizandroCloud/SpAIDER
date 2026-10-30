function [c, ceq] = mycon_reactor63_haar(w, m, Hi, iDWT, p)
% Restricoes - Exemplo 6.3 (Srinivasan, Palanki & Bonvin, 2003)
%
% w = [xdA(1:m); xdB(1:m); u(1:m)]
%
% Diferente do Exemplo 6.1: aqui as restricoes de concentracao sao
% TERMINAIS (so no final, nao em todo ponto de colocacao) --
% c so precisa checar o ULTIMO no, nao todos os m nos.
%
% Reconstrucao via Hi/iDWT (tf FIXO, ao contrario do 6.1):
%   cA(tau) = cAo + (xdA*inv(Hi))*iDWT
%   cB(tau) = cBo + (xdB*inv(Hi))*iDWT
%   V(tau)  = Vo  + tf*((u*inv(Hi))*iDWT)
%
% ceq: residuo da EDO para cA e cB (fator tf, pois d/dtau = tf*d/dt)
% c:   restricoes TERMINAIS (desigualdade, <=0):
%        cB(tf) <= cBf,max
%        cD(tf) <= cDf,max

    xdA = w(1     : m);
    xdB = w(m+1   : 2*m);
    u   = w(2*m+1 : 3*m);

    Hi_inv = inv(Hi); %#ok<MINV>

    cA = p.cAo + (xdA * Hi_inv) * iDWT;
    cB = p.cBo + (xdB * Hi_inv) * iDWT;
    V  = p.Vo  + p.tf * ((u * Hi_inv) * iDWT);

    dcA = -p.k1*cA.*cB - (u./V).*cA;
    dcB = -p.k1*cA.*cB - 2*p.k2*cB.^2 + (u./V).*(p.cBin - cB);

    ceq = [xdA - p.tf*dcA, xdB - p.tf*dcB];

    cA_f = cA(end); cB_f = cB(end); V_f = V(end);
    cD_f = ((cA_f + p.cBin - cB_f)*V_f - (p.cAo + p.cBin - p.cBo)*p.Vo) / (2*V_f);

    g_cB = cB_f - p.cBf_max;   % <=0
    g_cD = cD_f - p.cDf_max;   % <=0

    c = [g_cB, g_cD];

end
