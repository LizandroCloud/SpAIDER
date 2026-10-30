% Lizandro de Sousa Santos - Doutorado - PEQ - COPPE
% Exemplo 6.4 (Srinivasan, Palanki & Bonvin, 2003) - versao Haar (Hi/iDWT)
% "Non-isothermal semi-batch reactor with series reactions and a
%  heat removal constraint"
%
% Reacoes: A+B -> C -> D. Semi-batch, exotermico, NAO isotermico.
% Objetivo: maximizar nC(tf) = cC(tf)*V(tf).
% Variaveis manipuladas: DUAS -- vazao de B (u) e temperatura (T).
% Restricoes: limites de u e T; taxa de calor qrx<=qrx,max (CAMINHO);
% V(tf)<=Vmax (terminal). tf FIXO (0.5h).
%
% w = [xd1(1:m); xd2(1:m); u(1:m); T(1:m)]
%
% Requer no path: H_M.m, P_M.m, wav_mat.m, haar.m, haar2.m
%
% Parametros de temperatura/ΔH validados numericamente (cBo recalculado
% bate com a Tabela 6 do artigo assumindo ΔH1=-30000, ΔH2=-10000 J/mol
% e T em Kelvin = T[C]+273.15 dentro da Lei de Arrhenius).
%
% Referencia de otimo: nC*(tf) = 2.02 mol,
%   tempos de troca: tT = 0.05 h (T: Tmax->Tsens), tu = 0.3165 h (u: upath->umin)

clear; close all; clc;

%**************************************************
% Parametros (Tabela 6 do artigo)
%**************************************************
p.k1o = 4.0;      p.k2o = 800.0;
p.E1  = 6000.0;   p.E2  = 20000.0;
p.R   = 8.31;
p.DH1 = -30000.0; p.DH2 = -10000.0;   % exotermicas
p.umin = 0.0;     p.umax = 1.0;
p.Tmin = 20.0;    p.Tmax = 50.0;      % Celsius
p.Vmax = 1.1;
p.qrx_max = 1.5e5;
p.cAo = 10.0; p.cBo = 1.1685; p.cCo = 0.0; p.Vo = 1.0; p.cBin = 20.0;
p.tf  = 0.5;      % FIXO

p.x1o = p.Vo*p.cAo;             % = 10
p.x2o = p.Vo*(p.cAo + p.cCo);   % = 10

lambda_u = 0.3;
lambda_T = 0.0005;   % ver nota em myfun_reactor64_haar.m sobre sensibilidade de T

optNLP = optimset('Algorithm', 'sqp', 'Display', 'off', ...
    'TolX', 1e-11, 'TolFun', 1e-11, 'TolCon', 1e-11, ...
    'MaxFunEval', 200000, 'MaxIter', 6000);

%**************************************************
% Malha e chute estruturado (T alto no inicio, u alto no inicio)
%**************************************************
J = 5;
m = 2^(J+1);

[H, ~, Hi] = H_M(J);
[~, ~, PH] = P_M(J, H);
[~, iDWT]  = wav_mat(PH, Hi, J);

t_nodes0 = ((1:m) - 0.5) / m * p.tf;
u0 = zeros(1, m); T0 = zeros(1, m);
for i = 1:m
    t = t_nodes0(i);
    if t < 0.3165, u0(i) = 0.4; else, u0(i) = 0.0; end
    if t < 0.05,   T0(i) = p.Tmax; else, T0(i) = 30.0; end
end

w0 = [zeros(1,m), zeros(1,m), u0, T0];
lb = [-Inf(1,m), -Inf(1,m), p.umin*ones(1,m), p.Tmin*ones(1,m)];
ub = [ Inf(1,m),  Inf(1,m), p.umax*ones(1,m), p.Tmax*ones(1,m)];

[w, fval, exitflag, output] = fmincon(@(w) myfun_reactor64_haar(w, m, Hi, iDWT, p, lambda_u, lambda_T), w0, ...
    [], [], [], [], lb, ub, @(w) mycon_reactor64_haar(w, m, Hi, iDWT, p), optNLP);

Hi_inv = inv(Hi); %#ok<MINV>
x1 = p.x1o + (w(1:m)*Hi_inv)*iDWT;
x2 = p.x2o + (w(m+1:2*m)*Hi_inv)*iDWT;
x3 = p.Vo  + p.tf*((w(2*m+1:3*m)*Hi_inv)*iDWT);
nC_final = x2(end) - x1(end);

fprintf('J=%d (m=%d): nC(tf) = %.4f mol  | exitflag=%d | V(tf)=%.4f (Vmax=%.1f)\n', ...
    J, m, nC_final, exitflag, x3(end), p.Vmax);
fprintf('\nReferencia (Srinivasan et al., 2003): nC*(tf) = 2.02 mol\n');
fprintf('Tempos de troca de referencia: tT=0.05h (T), tu=0.3165h (u)\n');

%**************************************************
% Plot
%**************************************************
t_nodes = ((1:m) - 0.5) / m * p.tf;
u = w(2*m+1:3*m); T = w(3*m+1:4*m);

[c, ceq] = mycon_reactor64_haar(w, m, Hi, iDWT, p);
qrx = c(1:m) + p.qrx_max;   % desfaz o "-qrx_max" p/ recuperar qrx(t)

figure;
subplot(3,1,1);
stairs(t_nodes, u, 'LineWidth', 1.2);
ylabel('u [l/h]'); title(sprintf('Exemplo 6.4 (Haar) - nC(tf)=%.4f mol', nC_final));

subplot(3,1,2);
plot(t_nodes, T, '-', 'LineWidth', 1.2);
ylabel('T [C]');

subplot(3,1,3);
plot(t_nodes, qrx, '-', 'LineWidth', 1.2); hold on; yline(p.qrx_max, '--k');
xlabel('t [h]'); ylabel('q_{rx} [J/h]');
legend('q_{rx}(t)', 'q_{rx,max}', 'Location', 'best');
