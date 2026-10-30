function f = myfun_reactor63_haar(w, m, Hi, iDWT, p, lambda_reg)
% Funcao objetivo - Exemplo 6.3 (Srinivasan, Palanki & Bonvin, 2003)
% "Isothermal semi-batch reactor with parallel reactions and
%  selectivity constraints"
%
% w = [xdA(1:m); xdB(1:m); u(1:m)]
%   xdA, xdB = amostras de cA_dot, cB_dot (tempo normalizado tau=[0,1])
%   u        = vazao de B amostrada nos mesmos pontos
%
% tf e' FIXO aqui (p.tf = 250 min) -- diferente do Exemplo 6.1.
%
% max V(tf)*cC(tf) = nC(tf) = cAo*Vo - cA(tf)*V(tf)
% min f = -nC(tf) + lambda_reg*sum((u_k+1-u_k)^2)

    if nargin < 6
        lambda_reg = 5000;
        % OBS: a escala de u aqui e' ~1e-3 (umax=0.001), entao
        % diff(u)^2 e' ~1e-7 -- um lambda pequeno (ex: 0.01, que
        % funcionou no Exemplo 6.1 onde u~0.1) fica sem efeito aqui.
        % lambda=5000 foi calibrado numericamente para produzir um
        % perfil de u liso sem alterar nC(tf) (validado em Python).
    end

    xdA = w(1 : m);
    u   = w(2*m+1 : 3*m);

    Hi_inv = inv(Hi); %#ok<MINV>
    cA = p.cAo + (xdA * Hi_inv) * iDWT;
    V  = p.Vo  + p.tf * ((u * Hi_inv) * iDWT);

    nC_final = p.cAo*p.Vo - cA(end)*V(end);

    f = -nC_final + lambda_reg * sum(diff(u).^2);

end
