function [f, gradf] = myfun_reactor63_haar_grad(w, m, Hi, iDWT, p, lambda_reg)
% Funcao objetivo COM GRADIENTE ANALITICO - Exemplo 6.3, versao Haar
% w = [xdA(1:m); xdB(1:m); u(1:m)]   (tf e' FIXO aqui, p.tf, nao esta em w)
%
% min f = -nC(tf) + lambda_reg*sum((u_k+1-u_k)^2)
%   nC(tf) = cAo*Vo - cA(tf)*V(tf)
%
% Gradiente validado contra diferencas finitas em Python (erro ~1e-9).

    if nargin < 6
        lambda_reg = 5000;
    end

    xdA = w(1     : m);
    u   = w(2*m+1 : 3*m);

    A = Hi \ iDWT;

    cA = p.cAo + xdA*A;
    V  = p.Vo  + p.tf*(u*A);

    nC_final = p.cAo*p.Vo - cA(end)*V(end);
    f = -nC_final + lambda_reg*sum(diff(u).^2);

    if nargout > 1
        n = 3*m;
        gradf = zeros(n, 1);

        % d(f)/d(xdA_j) = V(end)*A(j,end)
        gradf(1:m) = V(end) * A(:, end);

        % d(f)/d(xdB_j) = 0  (ja inicializado com zeros)

        % d(f)/d(u_j) = cA(end)*tf*A(j,end) + lambda_reg*d(regularizacao)/d(u_j)
        g_reg = zeros(m, 1);
        g_reg(1) = -2*(u(2)-u(1));
        g_reg(m) = 2*(u(m)-u(m-1));
        g_reg(2:m-1) = 2*(u(2:m-1)-u(1:m-2)) - 2*(u(3:m)-u(2:m-1));

        gradf(2*m+1 : 3*m) = cA(end)*p.tf*A(:, end) + lambda_reg*g_reg;
    end

end
