function [c, ceq, dc, dceq] = mycon_case63_haar(w, m, Hi, iDWT, x0)
% Restricoes - Caso 6.3, Bell & Sargent (2000), usando a matriz de Haar.
%
% w = [xd1(1:m), xd2(1:m), u(1:m)]
%   xd1, xd2 = amostras de x1_dot, x2_dot nos pontos de colocacao
%   u        = controle amostrado nos mesmos pontos
%
% Reconstroi x1(t), x2(t) via Hi/iDWT (igual myfun_vo.m):
%   DX1 = xd1*inv(Hi);  X1 = DX1*iDWT + x0(1);
%   DX2 = xd2*inv(Hi);  X2 = DX2*iDWT + x0(2);
%
% E entao impoe que a taxa amostrada (xd1, xd2) seja CONSISTENTE com a
% dinamica bilinear avaliada em X1, X2, u -- isso e' o "defeito" da EDO
% nos pontos de colocacao, analogo ao ceq comentado no seu mycon.m
% original (agora efetivamente usado, com a EDO nao-linear do problema):
%   dx1/dt = u*(10*x2 - x1)
%   dx2/dt = u*(x1 - 10*x2) - (1-u)*x2

    xd1 = w(1:m);
    xd2 = w(m+1:2*m);
    u   = w(2*m+1:3*m);

    DX1 = xd1 * inv(Hi);              %#ok<MINV>
    X1  = DX1 * iDWT + x0(1);

    DX2 = xd2 * inv(Hi);              %#ok<MINV>
    X2  = DX2 * iDWT + x0(2);

    rhs1 = u .* (10.0*X2 - X1);
    rhs2 = u .* (X1 - 10.0*X2) - (1.0 - u) .* X2;

    ceq = [xd1 - rhs1, xd2 - rhs2];
    c = [];
    dc = [];
    dceq = [];

end
