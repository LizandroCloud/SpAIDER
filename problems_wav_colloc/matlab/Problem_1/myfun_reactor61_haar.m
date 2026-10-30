function [f, gradf] = myfun_reactor61_haar_grad(w, m, lambda_reg)
% Funcao objetivo COM GRADIENTE ANALITICO - Exemplo 6.1, versao Haar
% w = [tf; xdA(1:m); xdB(1:m); u(1:m)]
%
% min tf + lambda_reg*sum((u_k+1-u_k)^2)
%
% O termo de regularizacao (mesmo motivo do myfun_reactor61_haar.m
% original) evita o serrilhado numerico em u -- sem ele, o gradiente
% analitico converge mais rapido mas para o MESMO tipo de solucao
% ruidosa que a versao sem gradiente, pois a causa raiz (direcao
% quase-plana em u) nao tem relacao com como o gradiente e' calculado.

    if nargin < 3
        lambda_reg = 0.01;
    end

    n = length(w);
    u = w(2+2*m : 1+3*m);

    f = w(1) + lambda_reg*sum(diff(u).^2);

    if nargout > 1
        gradf = zeros(n, 1);
        gradf(1) = 1.0;

        % gradiente de lambda_reg*sum(diff(u)^2) em relacao a u
        % (validado contra diferencas finitas, erro ~1e-10)
        g_reg = zeros(m, 1);
        g_reg(1) = -2*(u(2)-u(1));
        g_reg(m) = 2*(u(m)-u(m-1));
        g_reg(2:m-1) = 2*(u(2:m-1)-u(1:m-2)) - 2*(u(3:m)-u(2:m-1));

        gradf(2+2*m : 1+3*m) = lambda_reg * g_reg;
    end

end
