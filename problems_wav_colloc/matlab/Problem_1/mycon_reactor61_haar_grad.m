function [c, ceq, gradc, gradceq] = mycon_reactor61_haar_grad(w, m, Hi, iDWT, p)
% Restricoes COM GRADIENTE ANALITICO - Exemplo 6.1, versao Haar.
%
% w = [tf; xdA(1:m); xdB(1:m); u(1:m)]
%
% Formulas do gradiente derivadas via regra da cadeia sobre a
% reconstrucao Hi/iDWT e validadas contra diferencas finitas em Python
% (erro maximo ~1e-10). Ver verify_gradients_61.py para a derivacao.
%
% Notacao: A = inv(Hi)*iDWT (m x m, CONSTANTE -- nao depende de w).
%   cA = cAo + xdA*A     => d(cA_k)/d(xdA_j) = A(j,k)
%   cB = cBo + xdB*A     => d(cB_k)/d(xdB_j) = A(j,k)
%   V  = Vo + tf*(u*A)   => d(V_k)/d(tf) = (u*A)_k ; d(V_k)/d(u_j) = tf*A(j,k)
%
% ceq_k = [xdA_k - tf*f1_k ; xdB_k - tf*f2_k], com
%   f1 = -k*cA*cB - (u/V)*cA
%   f2 = -k*cA*cB + (u/V)*(cBin-cB)
%
% c = [g_nC ; g_cB(1:m) ; g_V(1:m)]

    tf  = w(1);
    xdA = w(2 : 1+m);
    xdB = w(2+m : 1+2*m);
    u   = w(2+2*m : 1+3*m);

    A = Hi \ iDWT;    % = inv(Hi)*iDWT, mais estavel numericamente que inv() explicito

    cA = p.cAo + xdA*A;
    cB = p.cBo + xdB*A;
    V  = p.Vo  + tf*(u*A);

    f1 = -p.k*cA.*cB - (u./V).*cA;
    f2 = -p.k*cA.*cB + (u./V).*(p.cBin - cB);

    ceq = [xdA - tf*f1, xdB - tf*f2];

    nC_final = p.cAo*p.Vo - cA(end)*V(end);
    g_nC = p.nCdes - nC_final;
    g_cB = cB - p.cBmax;
    g_V  = V  - p.Vmax;
    c = [g_nC, g_cB, g_V];

    if nargout > 2

        n = 1 + 3*m;
        At = A';   % (A transposta, m x m) -- usada varias vezes abaixo

        % ---------- derivadas parciais de f1, f2 (elemento a elemento) ----------
        df1_dcA = -p.k*cB - u./V;
        df1_dcB = -p.k*cA;
        df1_dV  = (u./V.^2).*cA;
        df1_du  = -cA./V;

        df2_dcA = -p.k*cB;
        df2_dcB = -p.k*cA - u./V;
        df2_dV  = -(u./V.^2).*(p.cBin - cB);
        df2_du  = (p.cBin - cB)./V;

        dV_dtf = (u*A)';   % (m x 1): d(V_k)/d(tf)

        I = eye(m);

        % ---------- Jacobianos de ceq1, ceq2 (linha k = d(ceq_k)/d(.)) ----------
        J_ceq1_xdA = I - tf*bsxfun(@times, df1_dcA', At);
        J_ceq1_xdB =    - tf*bsxfun(@times, df1_dcB', At);
        J_ceq1_u   =    - tf*diag(df1_du) - tf^2*bsxfun(@times, df1_dV', At);
        d_ceq1_dtf = -f1' - tf*(df1_dV' .* dV_dtf);

        J_ceq2_xdA =    - tf*bsxfun(@times, df2_dcA', At);
        J_ceq2_xdB = I - tf*bsxfun(@times, df2_dcB', At);
        J_ceq2_u   =    - tf*diag(df2_du) - tf^2*bsxfun(@times, df2_dV', At);
        d_ceq2_dtf = -f2' - tf*(df2_dV' .* dV_dtf);

        gradceq = zeros(n, 2*m);
        gradceq(1, 1:m)             = d_ceq1_dtf';
        gradceq(2:1+m, 1:m)         = J_ceq1_xdA';
        gradceq(2+m:1+2*m, 1:m)     = J_ceq1_xdB';
        gradceq(2+2*m:1+3*m, 1:m)   = J_ceq1_u';

        gradceq(1, m+1:2*m)           = d_ceq2_dtf';
        gradceq(2:1+m, m+1:2*m)       = J_ceq2_xdA';
        gradceq(2+m:1+2*m, m+1:2*m)   = J_ceq2_xdB';
        gradceq(2+2*m:1+3*m, m+1:2*m) = J_ceq2_u';

        % ---------- gradiente de c = [g_nC, g_cB(1:m), g_V(1:m)] ----------
        gradc = zeros(n, 1 + 2*m);

        % g_nC (coluna 1)
        gradc(1, 1) = cA(end)*dV_dtf(end);
        gradc(2:1+m, 1) = V(end) * A(:, end);
        gradc(2+2*m:1+3*m, 1) = cA(end)*tf*A(:, end);

        % g_cB (colunas 2..m+1)
        gradc(2+m:1+2*m, 2:1+m) = A;

        % g_V (colunas m+2..2m+1)
        gradc(1, m+2:2*m+1) = dV_dtf';
        gradc(2+2*m:1+3*m, m+2:2*m+1) = tf*A;

    end

end
