% Lizandro de Sousa Santos - Doutorado - PEQ - COPPE
% Exemplo 6.1 (Srinivasan, Palanki & Bonvin, 2003) - versao Haar (Hi/iDWT)
% COM GRADIENTE ANALITICO (myfun_reactor61_haar_grad.m / mycon_reactor61_haar_grad.m)
%
% Mesma fisica/formulacao do main_reactor61_haar.m, mas agora informando
% o gradiente exato ao fmincon via 'SpecifyObjectiveGradient' e
% 'SpecifyConstraintGradient', em vez de deixar o fmincon estimar tudo
% por diferencas finitas.
%
% Validado em Python (scipy, mesma formula): erro do gradiente analitico
% vs diferencas finitas ~1e-10, e o numero de avaliacoes de funcao caiu
% de ~16500 (sem gradiente) para ~800 (com gradiente) no mesmo problema
% -- ~20x menos avaliacoes, mesmo resultado final (tf*).
%
% Referencia de otimo: tf* = 19.80 h, ts = 11.44 h

clear; close all; clc;

%**************************************************
% Parametros (Tabela 3 do artigo)
%**************************************************
p.k    = 0.0482;
p.DH   = 60000.0;
p.rho  = 900.0;
p.cp   = 4.2;
p.cBin = 2.0;
p.umin = 0.0;
p.umax = 0.1;
p.Tmax = 80.0;
p.Vmax = 1.0;
p.nCdes = 0.6;
p.cAo  = 2.0;
p.Vo   = 0.7;
T      = 70.0;

p.cBmax = p.rho*p.cp*(p.Tmax - T) / p.DH;
p.cBo   = p.cBmax;

fprintf('cBmax = %.4f mol/l\n', p.cBmax);

lambda_reg = 0.1;   % peso da regularizacao de suavidade em u (ver myfun_reactor61_haar_grad.m)

% *** aqui esta a diferenca-chave: informar o fmincon que myfun/mycon
% fornecem gradiente analitico, em vez de deixa-lo estimar por
% diferencas finitas ***
% NOTA: 'GradObj'/'GradConstr' sao os nomes de opcao do optimset em
% versoes mais antigas do MATLAB (ex: R2015). Em versoes mais recentes
% (R2016a+), o optimoptions usa 'SpecifyObjectiveGradient' e
% 'SpecifyConstraintGradient' (com true/false em vez de 'on'/'off').
% Se sua versao reclamar do nome oposto, troque para o outro par.
optNLP = optimset('Algorithm', 'sqp', 'Display', 'off', ...
    'TolX', 1e-10, 'TolFun', 1e-10, 'TolCon', 1e-10, ...
    'MaxFunEval', 200000, 'MaxIter', 5000, ...
    'GradObj', 'on', ...
    'GradConstr', 'on');

%**************************************************
% Malha e chute inicial
%**************************************************
J = 6;
m = 2^(J+1);

[H, ~, Hi] = H_M(J);
[~, ~, PH] = P_M(J, H);
[~, iDWT]  = wav_mat(PH, Hi, J);

w0 = [19.8, zeros(1,m), zeros(1,m), 0.05*ones(1,m)];
lb = [0.1,  -Inf(1,m), -Inf(1,m), p.umin*ones(1,m)];
ub = [100,   Inf(1,m),  Inf(1,m), p.umax*ones(1,m)];

tic;
[w, fval, exitflag, output] = fmincon(@(w) myfun_reactor61_haar_grad(w, m, lambda_reg), w0, ...
    [], [], [], [], lb, ub, @(w) mycon_reactor61_haar_grad(w, m, Hi, iDWT, p), optNLP);
tempo_com_grad = toc;

fprintf('\nCOM gradiente analitico:\n');
fprintf('  tf* = %.4f h  | exitflag=%d | funcCount=%d | tempo=%.3f s\n', ...
    w(1), exitflag, output.funcCount, tempo_com_grad);

%**************************************************
% Comparacao: mesma otimizacao SEM gradiente analitico (para referencia)
%**************************************************
optNLP_semgrad = optimset(optNLP, 'GradObj', 'off', ...
    'GradConstr', 'off');

tic;
[w2, fval2, exitflag2, output2] = fmincon(@(w) myfun_reactor61_haar_grad(w, m, lambda_reg), w0, ...
    [], [], [], [], lb, ub, @(w) mycon_reactor61_haar_grad(w, m, Hi, iDWT, p), optNLP_semgrad);
tempo_sem_grad = toc;

fprintf('\nSEM gradiente analitico (diferencas finitas):\n');
fprintf('  tf* = %.4f h  | exitflag=%d | funcCount=%d | tempo=%.3f s\n', ...
    w2(1), exitflag2, output2.funcCount, tempo_sem_grad);

fprintf('\nRazao de avaliacoes de funcao (sem/com gradiente): %.1fx\n', ...
    output2.funcCount / output.funcCount);

fprintf('\nReferencia (Srinivasan et al., 2003): tf* = 19.80 h\n');

%**************************************************
% Plot da solucao (com gradiente)
%**************************************************
tf  = w(1);
xdA = w(2 : 1+m);
xdB = w(2+m : 1+2*m);
u   = w(2+2*m : 1+3*m);

A = Hi \ iDWT;
cA = p.cAo + xdA*A;
cB = p.cBo + xdB*A;
V  = p.Vo  + tf*(u*A);
cC = (p.cAo*p.Vo - cA.*V) ./ V;

t_nodes = ((1:m) - 0.5) / m * tf;

figure;
subplot(2,1,1);
stairs(t_nodes, u, 'LineWidth', 1.2);
ylabel('u [l/h]'); title(sprintf('Exemplo 6.1 (Haar, c/ gradiente analitico) - tf=%.2f h', tf));

subplot(2,1,2);
plot(t_nodes, cA, '-', t_nodes, cB, '-', t_nodes, cC, '-', 'LineWidth', 1.2);
hold on; yline(p.cBmax, '--k');
xlabel('t [h]'); ylabel('concentracao [mol/l]');
legend('c_A', 'c_B', 'c_C', 'c_{B,max}', 'Location', 'best');
