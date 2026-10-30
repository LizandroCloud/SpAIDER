function f = myfun_reactor64_haar(w, m, Hi, iDWT, p, lambda_u, lambda_T)
% Funcao objetivo - Exemplo 6.4 (Srinivasan, Palanki & Bonvin, 2003)
% "Non-isothermal semi-batch reactor with series reactions and a
%  heat removal constraint"
%
% w = [xd1(1:m); xd2(1:m); u(1:m); T(1:m)]
%   xd1, xd2 = amostras de x1_dot, x2_dot (tempo normalizado tau=[0,1])
%     x1 = V*cA,  x2 = V*(cA+cC)   (modelo reduzido, Eqs. 76-79)
%   u = vazao de B [l/h],  T = temperatura do reator [C]
%
% tf e' FIXO (p.tf = 0.5 h). DUAS variaveis manipuladas (u, T).
%
% max nC(tf) = cC(tf)*V(tf) = x2(tf) - x1(tf)
% min f = -nC(tf) + lambda_u*sum((u_k+1-u_k)^2) + lambda_T*sum((T_k+1-T_k)^2)
%
% NOTA: assim como no arco usens do Exemplo 6.3, o perfil de T aqui tem
% uma direcao com pouca sensibilidade de custo -- varios perfis de T
% dao nC(tf) parecido. lambda_T pequeno preserva melhor a forma
% qualitativa esperada (T alto no inicio, cai, sobe no final);
% lambda_T maior da um perfil mais liso mas mais achatado (validado
% numericamente).

    if nargin < 6, lambda_u = 0.3; end
    if nargin < 7, lambda_T = 0.0005; end

    xd1 = w(1 : m);
    xd2 = w(m+1 : 2*m);
    u   = w(2*m+1 : 3*m);
    T   = w(3*m+1 : 4*m);

    Hi_inv = inv(Hi); %#ok<MINV>
    x1 = p.x1o + (xd1 * Hi_inv) * iDWT;
    x2 = p.x2o + (xd2 * Hi_inv) * iDWT;

    nC_final = x2(end) - x1(end);

    f = -nC_final + lambda_u*sum(diff(u).^2) + lambda_T*sum(diff(T).^2);

end
