% Lizandro de Sousa Santos - Doutorado - PEQ - COPPE
% Exemplo 6.3 (Srinivasan, Palanki & Bonvin, 2003) - versao Haar (Hi/iDWT)
% "Isothermal semi-batch reactor with parallel reactions and
%  selectivity constraints"
%
% Reacoes: A+B -> C, 2B -> D. Semi-batch, isotermico.
% Objetivo: maximizar V(tf)*cC(tf) = nC(tf).
% Variavel manipulada: vazao de alimentacao de B.
% Restricoes: limites de u; cB(tf)<=cBf,max; cD(tf)<=cDf,max (SO no
% final, diferente do Exemplo 6.1 que tinha restricao de CAMINHO).
% tf e' FIXO aqui (250 min) -- diferente do Exemplo 6.1 (tf livre).
%
% w = [xdA(1:m); xdB(1:m); u(1:m)]
%   xdA,xdB  = amostras de cA_dot, cB_dot (tempo normalizado tau=[0,1])
%   u        = vazao de B amostrada nos mesmos pontos de colocacao
%
% Requer no path: H_M.m, P_M.m, wav_mat.m, haar.m, haar2.m
%
% Referencia de otimo: nC*(tf) = 0.43 mol,
%   tempos de troca: tm = 20.25 min (umax->usens), ts = 205 min (usens->umin)

clear; close all; clc;

%**************************************************
% Parametros (Tabela 5 do artigo)
%**************************************************
p.k1  = 0.053;
p.k2  = 0.128;
p.cBin = 5.0;
p.umin = 0.0;
p.umax = 0.001;
p.cBf_max = 0.025;
p.cDf_max = 0.15;
p.cAo = 0.72;
p.cBo = 0.05;
p.Vo  = 1.0;
p.tf  = 250.0;   % FIXO

lambda_reg = 5000;   % peso da regularizacao de suavidade em u (ver nota em myfun_reactor63_haar.m)

optNLP = optimset('Algorithm', 'sqp', 'Display', 'off', ...
    'TolX', 1e-6, 'TolFun', 1e-6, 'TolCon', 1e-6, ...
    'MaxFunEval', 200000, 'MaxIter', 5000);

%**************************************************
% Passo 1: malha grosseira (J=3 -> m=16)
%**************************************************
J = 3;
m = 2^(J+1);

[H, ~, Hi] = H_M(J);
[~, ~, PH] = P_M(J, H);
[~, iDWT]  = wav_mat(PH, Hi, J);

w0 = [zeros(1,m), zeros(1,m), 0.5*p.umax*ones(1,m)];
lb = [-Inf(1,m), -Inf(1,m), p.umin*ones(1,m)];
ub = [ Inf(1,m),  Inf(1,m), p.umax*ones(1,m)];

[w, fval, exitflag, output] = fmincon(@(w) myfun_reactor63_haar(w, m, Hi, iDWT, p, lambda_reg), w0, ...
    [], [], [], [], lb, ub, @(w) mycon_reactor63_haar(w, m, Hi, iDWT, p), optNLP);

xdA = w(1:m); V = p.Vo + p.tf*((w(2*m+1:3*m)*inv(Hi))*iDWT); %#ok<MINV>
cA = p.cAo + (xdA*inv(Hi))*iDWT; %#ok<MINV>
nC = p.cAo*p.Vo - cA(end)*V(end);
fprintf('J=%d (m=%d): nC(tf) = %.4f mol  | exitflag=%d\n', J, m, nC, exitflag);

%**************************************************
% Passo 2: refinamento de malha, reaproveitando a solucao anterior
%**************************************************
w_prev = w;
m_prev = m;
for J2 = 4:1:6
    m2 = 2^(J2+1);
    rep = m2 / m_prev;

    xdA_prev = w_prev(1 : m_prev);
    xdB_prev = w_prev(m_prev+1 : 2*m_prev);
    u_prev   = w_prev(2*m_prev+1 : 3*m_prev);

    w0 = [repelem(xdA_prev, rep), repelem(xdB_prev, rep), repelem(u_prev, rep)];

    [H2, ~, Hi2] = H_M(J2);
    [~, ~, PH2]  = P_M(J2, H2);
    [~, iDWT2]   = wav_mat(PH2, Hi2, J2);

    lb2 = [-Inf(1,m2), -Inf(1,m2), p.umin*ones(1,m2)];
    ub2 = [ Inf(1,m2),  Inf(1,m2), p.umax*ones(1,m2)];

    [w, fval, exitflag, output] = fmincon(@(w) myfun_reactor63_haar(w, m2, Hi2, iDWT2, p, lambda_reg), w0, ...
        [], [], [], [], lb2, ub2, @(w) mycon_reactor63_haar(w, m2, Hi2, iDWT2, p), optNLP);

    [c, ceq] = mycon_reactor63_haar(w, m2, Hi2, iDWT2, p);
    xdA2 = w(1:m2); u2 = w(2*m2+1:3*m2);
    Hi2_inv = inv(Hi2); %#ok<MINV>
    cA2 = p.cAo + (xdA2*Hi2_inv)*iDWT2;
    V2  = p.Vo + p.tf*((u2*Hi2_inv)*iDWT2);
    nC2 = p.cAo*p.Vo - cA2(end)*V2(end);

    fprintf('J=%d (m=%d): nC(tf) = %.4f mol | exitflag=%d | max|ceq|=%.2e | max(c)=%.2e | max|du|=%.5f\n', ...
        J2, m2, nC2, exitflag, max(abs(ceq)), max(c), max(abs(diff(u2))));

    w_prev = w;
    m_prev = m2;
end

fprintf('\nReferencia (Srinivasan et al., 2003): nC*(tf) = 0.43 mol\n');
fprintf('Tempos de troca de referencia: tm = 20.25 min, ts = 205 min\n');

%**************************************************
% Reconstrucao final e plot
%**************************************************
m = m_prev;
[H, ~, Hi] = H_M(log2(m) - 1);
[~, ~, PH] = P_M(log2(m) - 1, H);
[~, iDWT]  = wav_mat(PH, Hi, log2(m) - 1);

xdA = w(1:m); xdB = w(m+1:2*m); u = w(2*m+1:3*m);
Hi_inv = inv(Hi); %#ok<MINV>
cA = p.cAo + (xdA*Hi_inv)*iDWT;
cB = p.cBo + (xdB*Hi_inv)*iDWT;
V  = p.Vo  + p.tf*((u*Hi_inv)*iDWT);
cC = (p.cAo*p.Vo - cA.*V) ./ V;
cD = ((cA + p.cBin - cB).*V - (p.cAo + p.cBin - p.cBo)*p.Vo) ./ (2*V);

%**************************************************
% Comparacao com a expressao analitica u_sens (arco do meio, ver
% usens_analitico_63.m -- re-derivada simbolicamente, ja que a Eq. 71
% do artigo original saiu corrompida no OCR)
%**************************************************
u_sens_analitico = usens_analitico_63(cA, cB, V, p.k1, p.k2, p.cBin);

mid = (u > 1e-6) & (u < p.umax - 1e-6);   % arco do meio (nem umax nem umin)
if any(mid)
    erro_usens = max(abs(u(mid) - u_sens_analitico(mid)));
    fprintf('max|u_numerico - u_sens_analitico| no arco do meio = %.3e\n', erro_usens);
end

t_nodes = ((1:m) - 0.5) / m * p.tf;

figure;
subplot(2,1,1);
stairs(t_nodes, u, 'LineWidth', 1.2); hold on;
plot(t_nodes(mid), u_sens_analitico(mid), 'r--', 'LineWidth', 1.0);
ylabel('u [l/min]'); title(sprintf('Exemplo 6.3 (Haar) - vazao de B (nC(tf)=%.4f mol)', p.cAo*p.Vo - cA(end)*V(end)));
legend('u numerico (fmincon)', 'u_{sens} analitico', 'Location', 'best');

subplot(2,1,2);
plot(t_nodes, cB, '-', t_nodes, cD, '-', 'LineWidth', 1.2);
hold on; yline(p.cBf_max, '--k'); yline(p.cDf_max, '--r');
xlabel('t [min]'); ylabel('concentracao [mol/l]');
legend('c_B', 'c_D', 'c_{B,max}', 'c_{D,max}', 'Location', 'best');
