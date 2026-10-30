% Lizandro de Sousa Santos - Doutorado - PEQ - COPPE
% Exemplo 6.3 (Srinivasan, Palanki & Bonvin, 2003) - versao Haar (Hi/iDWT)
% COM GRADIENTE ANALITICO (myfun_reactor63_haar_grad.m / mycon_reactor63_haar_grad.m)
%
% Mesma formulacao do main_reactor63_haar.m, mas informando o gradiente
% exato ao fmincon (GradObj/GradConstr 'on' -- sintaxe do optimset
% compativel com MATLAB anteriores ao R2016a; troque para
% 'SpecifyObjectiveGradient'/'SpecifyConstraintGradient' com
% true/false se sua versao for mais recente e reclamar do nome).
%
% Gradiente validado em Python contra diferencas finitas (erro ~1e-9).
%
% Referencia de otimo: nC*(tf) = 0.43 mol

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

lambda_reg = 5000;   % ver nota em myfun_reactor63_haar_grad.m sobre escala de u

% *** GradObj/GradConstr 'on' informa ao fmincon que myfun/mycon
% fornecem gradiente analitico (2o/3o-4o argumentos de saida) ***
optNLP = optimset('Algorithm', 'sqp', 'Display', 'off', ...
    'TolX', 1e-12, 'TolFun', 1e-12, 'TolCon', 1e-12, ...
    'MaxFunEval', 200000, 'MaxIter', 5000, ...
    'GradObj', 'on', 'GradConstr', 'on');

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

tic;
[w, fval, exitflag, output] = fmincon(@(w) myfun_reactor63_haar_grad(w, m, Hi, iDWT, p, lambda_reg), w0, ...
    [], [], [], [], lb, ub, @(w) mycon_reactor63_haar_grad(w, m, Hi, iDWT, p), optNLP);
tempo_com_grad = toc;

A = Hi \ iDWT;
cA = p.cAo + w(1:m)*A;
V  = p.Vo  + p.tf*(w(2*m+1:3*m)*A);
nC = p.cAo*p.Vo - cA(end)*V(end);
fprintf('J=%d (m=%d): nC(tf) = %.4f mol  | exitflag=%d | funcCount=%d | tempo=%.3f s\n', ...
    J, m, nC, exitflag, output.funcCount, tempo_com_grad);

%**************************************************
% Passo 2: refinamento de malha, reaproveitando a solucao anterior
%**************************************************
w_prev = w;
m_prev = m;
for J2 = [4 : 1: 6]
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

    [w, fval, exitflag, output] = fmincon(@(w) myfun_reactor63_haar_grad(w, m2, Hi2, iDWT2, p, lambda_reg), w0, ...
        [], [], [], [], lb2, ub2, @(w) mycon_reactor63_haar_grad(w, m2, Hi2, iDWT2, p), optNLP);

    A2 = Hi2 \ iDWT2;
    xdA2 = w(1:m2); u2 = w(2*m2+1:3*m2);
    cA2 = p.cAo + xdA2*A2;
    V2  = p.Vo + p.tf*(u2*A2);
    nC2 = p.cAo*p.Vo - cA2(end)*V2(end);

    [c, ceq] = mycon_reactor63_haar_grad(w, m2, Hi2, iDWT2, p);
    fprintf('J=%d (m=%d): nC(tf) = %.4f mol | exitflag=%d | funcCount=%d | max|ceq|=%.2e | max(c)=%.2e | max|du|=%.5f\n', ...
        J2, m2, nC2, exitflag, output.funcCount, max(abs(ceq)), max(c), max(abs(diff(u2))));

    w_prev = w;
    m_prev = m2;
end

fprintf('\nReferencia (Srinivasan et al., 2003): nC*(tf) = 0.43 mol\n');

%**************************************************
% Comparacao: mesma otimizacao (malha final) SEM gradiente analitico
%**************************************************
optNLP_semgrad = optimset(optNLP, 'GradObj', 'off', 'GradConstr', 'off');

m2 = m_prev;
rep = m2 / m; % chute inicial da malha grosseira, igual ao passo 1 original
w0 = [zeros(1,m2), zeros(1,m2), 0.5*p.umax*ones(1,m2)];
lb2 = [-Inf(1,m2), -Inf(1,m2), p.umin*ones(1,m2)];
ub2 = [ Inf(1,m2),  Inf(1,m2), p.umax*ones(1,m2)];
[H2, ~, Hi2] = H_M(log2(m2)-1);
[~, ~, PH2]  = P_M(log2(m2)-1, H2);
[~, iDWT2]   = wav_mat(PH2, Hi2, log2(m2)-1);

tic;
[w_ng, ~, exitflag_ng, output_ng] = fmincon(@(w) myfun_reactor63_haar_grad(w, m2, Hi2, iDWT2, p, lambda_reg), w0, ...
    [], [], [], [], lb2, ub2, @(w) mycon_reactor63_haar_grad(w, m2, Hi2, iDWT2, p), optNLP_semgrad);
tempo_sem_grad = toc;

fprintf('\nSEM gradiente analitico (mesma malha final, m=%d):\n', m2);
fprintf('  exitflag=%d | funcCount=%d | tempo=%.3f s\n', exitflag_ng, output_ng.funcCount, tempo_sem_grad);
fprintf('Razao de avaliacoes de funcao (sem/com gradiente): %.1fx\n', output_ng.funcCount / output.funcCount);

%**************************************************
% Plot da solucao final (com gradiente)
%**************************************************
m = m_prev;
A = Hi2 \ iDWT2;
xdA = w(1:m); xdB = w(m+1:2*m); u = w(2*m+1:3*m);
cA = p.cAo + xdA*A;
cB = p.cBo + xdB*A;
V  = p.Vo  + p.tf*(u*A);
cD = ((cA + p.cBin - cB).*V - (p.cAo + p.cBin - p.cBo)*p.Vo) ./ (2*V);

t_nodes = ((1:m) - 0.5) / m * p.tf;

figure;
subplot(2,1,1);
stairs(t_nodes, u, 'LineWidth', 1.2);
ylabel('u [l/min]'); title(sprintf('Exemplo 6.3 (Haar, c/ gradiente analitico) - nC(tf)=%.4f mol', p.cAo*p.Vo - cA(end)*V(end)));

subplot(2,1,2);
plot(t_nodes, cB, '-', t_nodes, cD, '-', 'LineWidth', 1.2);
hold on; yline(p.cBf_max, '--k'); yline(p.cDf_max, '--r');
xlabel('t [min]'); ylabel('concentracao [mol/l]');
legend('c_B', 'c_D', 'c_{B,max}', 'c_{D,max}', 'Location', 'best');
