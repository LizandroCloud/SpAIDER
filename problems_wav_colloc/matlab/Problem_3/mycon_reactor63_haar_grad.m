function [c, ceq, gradc, gradceq] = mycon_reactor63_haar_grad(w, m, Hi, iDWT, p)
% Restricoes COM GRADIENTE ANALITICO - Exemplo 6.3, versao Haar.
%
% w = [xdA(1:m); xdB(1:m); u(1:m)]   (tf FIXO, p.tf)
%
% ceq_k = [xdA_k - tf*f1_k ; xdB_k - tf*f2_k], com
%   f1 = -k1*cA*cB - (u/V)*cA
%   f2 = -k1*cA*cB - 2*k2*cB^2 + (u/V)*(cBin-cB)
%
% c = [g_cB, g_cD]  (restricoes TERMINAIS, so no ultimo no)
%
% Gradiente validado contra diferencas finitas em Python (erro ~1e-9).
% Usa bsxfun (nao "implicit expansion") para compatibilidade com
% MATLAB anteriores ao R2016b.

    xdA = w(1     : m);
    xdB = w(m+1   : 2*m);
    u   = w(2*m+1 : 3*m);

    A = Hi \ iDWT;

    cA = p.cAo + xdA*A;
    cB = p.cBo + xdB*A;
    V  = p.Vo  + p.tf*(u*A);

    f1 = -p.k1*cA.*cB - (u./V).*cA;
    f2 = -p.k1*cA.*cB - 2*p.k2*cB.^2 + (u./V).*(p.cBin - cB);

    ceq = [xdA - p.tf*f1, xdB - p.tf*f2];

    cA_f = cA(end); cB_f = cB(end); V_f = V(end);
    cD_f = ((cA_f + p.cBin - cB_f)*V_f - (p.cAo + p.cBin - p.cBo)*p.Vo) / (2*V_f);

    g_cB = cB_f - p.cBf_max;
    g_cD = cD_f - p.cDf_max;
    c = [g_cB, g_cD];

    if nargout > 2

        n = 3*m;
        At = A';

        df1_dcA = -p.k1*cB - u./V;
        df1_dcB = -p.k1*cA;
        df1_dV  = (u./V.^2).*cA;
        df1_du  = -cA./V;

        df2_dcA = -p.k1*cB;
        df2_dcB = -p.k1*cA - 4*p.k2*cB - u./V;
        df2_dV  = -(u./V.^2).*(p.cBin - cB);
        df2_du  = (p.cBin - cB)./V;

        I = eye(m);

        J_ceq1_xdA = I - p.tf*bsxfun(@times, df1_dcA', At);
        J_ceq1_xdB =    - p.tf*bsxfun(@times, df1_dcB', At);
        J_ceq1_u   =    - p.tf*diag(df1_du) - p.tf^2*bsxfun(@times, df1_dV', At);

        J_ceq2_xdA =    - p.tf*bsxfun(@times, df2_dcA', At);
        J_ceq2_xdB = I - p.tf*bsxfun(@times, df2_dcB', At);
        J_ceq2_u   =    - p.tf*diag(df2_du) - p.tf^2*bsxfun(@times, df2_dV', At);

        gradceq = zeros(n, 2*m);
        gradceq(1:m, 1:m)         = J_ceq1_xdA';
        gradceq(m+1:2*m, 1:m)     = J_ceq1_xdB';
        gradceq(2*m+1:3*m, 1:m)   = J_ceq1_u';
        gradceq(1:m, m+1:2*m)     = J_ceq2_xdA';
        gradceq(m+1:2*m, m+1:2*m) = J_ceq2_xdB';
        gradceq(2*m+1:3*m, m+1:2*m) = J_ceq2_u';

        % ---- gradiente de c = [g_cB, g_cD] ----
        gradc = zeros(n, 2);

        gradc(m+1:2*m, 1) = A(:, end);   % g_cB

        Kc = (p.cAo + p.cBin - p.cBo) * p.Vo;
        gradc(1:m, 2)       =  0.5 * A(:, end);
        gradc(m+1:2*m, 2)   = -0.5 * A(:, end);
        gradc(2*m+1:3*m, 2) = (Kc / (2*V_f^2)) * p.tf * A(:, end);

    end

end
