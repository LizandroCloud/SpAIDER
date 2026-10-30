function [c, ceq] = mycon_reactor64_haar(w, m, Hi, iDWT, p)
% Restricoes - Exemplo 6.4 (Srinivasan, Palanki & Bonvin, 2003)
%
% w = [xd1(1:m); xd2(1:m); u(1:m); T(1:m)]
%
% Modelo reduzido (Eqs. 76-79, com a correcao de sinal x2=V*(cA+cC) --
% ver nota abaixo):
%   x1_dot = -k1*x1*cB
%   x2_dot = k2*(x1-x2)
%   x3_dot = u                          (V, reconstruido direto)
%   cB = cBin + x1/x3 + (Vo/x3)*(cBo-cBin-cAo)
%   cC = (x2-x1)/x3
%   k1 = k1o*exp(-E1/(R*(T+273.15)))    (T em Celsius no vetor de decisao,
%                                         convertido para Kelvin aqui)
%   k2 = k2o*exp(-E2/(R*(T+273.15)))
%
% NOTA IMPORTANTE: a Eq. (76)-(79) do artigo, no PDF fornecido, saiu com
% "x2=V(cA/cC)" no OCR -- ambiguo entre (cA-cC) e (cA+cC). Verificou-se
% por derivacao direta (d/dt de V*cA e V*cC) que so x2=V*(cA+cC) produz
% a EDO x2_dot=k2*(x1-x2) declarada no artigo (a escolha (cA-cC) gera um
% termo espurio em k1). Usa-se aqui x2=V*(cA+cC), cC=(x2-x1)/x3.
%
% qrx = (-DH1)*k1*x1*cB + (-DH2)*k2*(x2-x1)   (o V cancela algebricamente
%        de (-DH1)*k1*cA*cB*V + (-DH2)*k2*cC*V, ver Eq. 80)
%
% ceq: residuo da EDO para x1, x2
% c:   qrx(t) <= qrx_max  (CAMINHO, todos os nos)  e  V(tf) <= Vmax (terminal)

    xd1 = w(1     : m);
    xd2 = w(m+1   : 2*m);
    u   = w(2*m+1 : 3*m);
    T   = w(3*m+1 : 4*m);

    Hi_inv = inv(Hi); %#ok<MINV>

    x1 = p.x1o + (xd1 * Hi_inv) * iDWT;
    x2 = p.x2o + (xd2 * Hi_inv) * iDWT;
    x3 = p.Vo  + p.tf * ((u * Hi_inv) * iDWT);

    TK = T + 273.15;
    k1 = p.k1o * exp(-p.E1 ./ (p.R*TK));
    k2 = p.k2o * exp(-p.E2 ./ (p.R*TK));
    cB = p.cBin + x1./x3 + (p.Vo./x3)*(p.cBo - p.cBin - p.cAo);

    ceq = [xd1 - p.tf*(-k1.*x1.*cB), xd2 - p.tf*(k2.*(x1-x2))];

    qrx = (-p.DH1)*k1.*x1.*cB + (-p.DH2)*k2.*(x2-x1);
    g_qrx = qrx - p.qrx_max;      % <=0, todos os nos
    g_V   = x3(end) - p.Vmax;     % <=0, terminal

    c = [g_qrx, g_V];

end
