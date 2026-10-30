function f = myfun_case63_haar(w, m, Hi, iDWT, x0)
% Funcao objetivo - Caso 6.3, Bell & Sargent (2000)
% Usa a matriz de Haar (Hi) e o operador de integracao (iDWT), no MESMO
% papel que exercem no seu myfun_vo.m: xd1/xd2 (as duas primeiras partes
% de w) sao AMOSTRAS de x1_dot e x2_dot nos pontos de colocacao de Haar
% -- nao coeficientes, nao valores de x1/x2 diretamente. Hi/iDWT fazem a
% conversao amostra->coeficiente->integral->amostra, exatamente como em
%   DX1 = xd*inv(Hi);  X1 = DX1*iDWT + x0;

    xd1 = w(1:m);
    xd2 = w(m+1:2*m);

    DX1 = xd1 * inv(Hi);              %#ok<MINV>
    X1  = DX1 * iDWT + x0(1);

    DX2 = xd2 * inv(Hi);              %#ok<MINV>
    X2  = DX2 * iDWT + x0(2);

    x3f = 1.0 - X1(end) - X2(end);
    f = -x3f;

end
