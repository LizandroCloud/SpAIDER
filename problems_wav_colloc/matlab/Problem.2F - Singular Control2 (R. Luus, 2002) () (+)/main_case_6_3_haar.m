% Lizandro de Sousa Santos - Doutorado - PEQ - COPPE
% Caso 6.3 (Bell & Sargent, 2000) - "mixed catalyst problem"
% Versao usando de fato a matriz de Haar (H_M.m, P_M.m, wav_mat.m) --
% mesma infraestrutura do seu main_vo.m / myfun_vo.m / mycon.m.
%
% Decisao w = [xd1(1:m), xd2(1:m), u(1:m)]:
%   xd1, xd2 = amostras de x1_dot, x2_dot nos pontos de colocacao de Haar
%   u        = controle amostrado nos mesmos pontos
% x1(t), x2(t) sao reconstruidos dentro de myfun_case63_haar.m /
% mycon_case63_haar.m via Hi/iDWT (nao sao variaveis de decisao
% diretas -- isso e' o que muda em relacao a versao "dominio fisico"
% que te mandei antes).
%
% Requer no path: H_M.m, P_M.m, wav_mat.m, haar.m, haar2.m (seus
% arquivos originais).
%
% Sistema:
%   dx1/dt = u*(10*x2 - x1)
%   dx2/dt = u*(x1 - 10*x2) - (1-u)*x2
%   x3(t)  = 1 - x1(t) - x2(t)
%   x1(0)=1.0, x2(0)=0.0,  0<=u(t)<=1,  t em [0,1]
%   max x3(1)  <=>  min f = -x3(1)
%
% Referencias de otimo:
%   Bell & Sargent (2000):  x3(1) = 0.048080
%   Vassiliadis (1993):     x3(1) = 0.0480557

clear; close all; clc;

x0 = [1.0, 0.0];

optNLP = optimset('Algorithm', 'sqp', 'Display', 'off', ...
    'TolX', 1e-10, 'TolFun', 1e-10, 'TolCon', 1e-10, 'MaxFunEval', 100000);

%**************************************************
% Passo 1: multi-start numa malha grosseira (J=2 -> m=8)
% (arco singular -> chute inicial importa; usamos varios)
%**************************************************
J = 2;
m = 2^(J+1);

[H, ~, Hi] = H_M(J);
[~, ~, PH] = P_M(J, H);
[~, iDWT]  = wav_mat(PH, Hi, J);

lb = [-Inf(1, 2*m), zeros(1, m)];
ub = [ Inf(1, 2*m), ones(1, m)];

best_fval = Inf;
best_w = [];
n_starts = 8;
rng(0);
for trial = 1:n_starts
    if trial == 1
        u0 = 0.17 * ones(1, m);
    else
        u0 = rand(1, m);
    end
    w0 = [zeros(1, m), zeros(1, m), u0];   % xd1=xd2=0 como chute (estado ~ constante)

    [w, fval] = fmincon(@(w) myfun_case63_haar(w, m, Hi, iDWT, x0), w0, ...
        [], [], [], [], lb, ub, ...
        @(w) mycon_case63_haar(w, m, Hi, iDWT, x0), optNLP);

    if fval < best_fval
        best_fval = fval;
        best_w = w;
    end
end
fprintf('J=%d (m=%d): x3(1) = %.6f\n', J, m, -best_fval);

%**************************************************
% Passo 2: refinamento de malha (J crescente), reaproveitando a solucao
%**************************************************
w_prev = best_w;
J_prev = J;
m_prev = m;
for J2 = [3 4 5]
    m2 = 2^(J2+1);
    rep = m2 / m_prev;

    xd1p = w_prev(1:m_prev);
    xd2p = w_prev(m_prev+1:2*m_prev);
    up   = w_prev(2*m_prev+1:3*m_prev);

    w0 = [repelem(xd1p, rep), repelem(xd2p, rep), repelem(up, rep)];

    [H2, ~, Hi2] = H_M(J2);
    [~, ~, PH2]  = P_M(J2, H2);
    [~, iDWT2]   = wav_mat(PH2, Hi2, J2);

    lb2 = [-Inf(1, 2*m2), zeros(1, m2)];
    ub2 = [ Inf(1, 2*m2), ones(1, m2)];

    [w, fval] = fmincon(@(w) myfun_case63_haar(w, m2, Hi2, iDWT2, x0), w0, ...
        [], [], [], [], lb2, ub2, ...
        @(w) mycon_case63_haar(w, m2, Hi2, iDWT2, x0), optNLP);

    fprintf('J=%d (m=%d): x3(1) = %.6f\n', J2, m2, -fval);

    w_prev = w;
    m_prev = m2;
    J_prev = J2; %#ok<NASGU>
end

fprintf('\nReferencia (Bell & Sargent, 2000): 0.048080\n');
fprintf('Referencia (Vassiliadis, 1993):    0.0480557\n');

%**************************************************
% Plot da solucao final: reconstroi x1(t), x2(t) via Hi/iDWT
%**************************************************
m = m_prev;
[H, ~, Hi] = H_M(log2(m) - 1);
[~, ~, PH] = P_M(log2(m) - 1, H);
[~, iDWT]  = wav_mat(PH, Hi, log2(m) - 1);

xd1 = w_prev(1:m);
xd2 = w_prev(m+1:2*m);
uf  = w_prev(2*m+1:3*m);

X1 = (xd1 * inv(Hi)) * iDWT + x0(1); %#ok<MINV>
X2 = (xd2 * inv(Hi)) * iDWT + x0(2); %#ok<MINV>

t_nodes = ((1:m) - 0.5) / m;   % pontos de colocacao de Haar

figure;
subplot(2,1,1);
stairs(t_nodes, uf, 'LineWidth', 1.2);
ylabel('u(t)'); title(sprintf('Caso 6.3 - colocacao via Haar (Hi/iDWT), m=%d', m));
ylim([-0.05 1.05]);

subplot(2,1,2);
plot(t_nodes, X1, '-o', 'MarkerSize', 3); hold on;
plot(t_nodes, X2, '-o', 'MarkerSize', 3);
xlabel('t'); ylabel('estados'); legend('x1', 'x2', 'Location', 'best');
