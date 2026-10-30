%RUN MODELO KADAM, 2004
clear all; clc;
%===================
% DADOS ENTRADA
Fenz = 0.0001;      % Vazão de enzima (L/min)
Fsol = 0.0004;      % Vazão de sólidos (L/min)
Fout = 0.0;         % Vazão de saida do reator (L/min)
Eenz = 1.937*(5/1); % [] enzima entra mixer (g/L)
Csol = 80*(5/4);    % [] celulose entra mixer (g/L)
Gsol = 0*(5/4);     % [] glicose entra mixer (g/L)
Fin = Fsol + Fenz;
Cin = (Fsol/Fin)*Csol;
Gin = (Fsol/Fin)*Gsol;
Ein = (Fenz/Fin)*Eenz;

%C.I. para EDO's
Cini = 80;      %(g/L)
Gini = 0;       %(g/L)
Eini = 1.937;   %(g/L)
Vini = 0.1;     %(L)
%===================

Mini= [Cini; Gini; Eini; Vini];
%Mini = [Cini; Gini];
tini = 0;
tfin = 196*60; %(min)
op = odeset('RelTol',1E-6, 'AbsTol',1E-8, 'BDF', 'on');

%Solve EDO's
[t,M] = ode15s('Avicel_mod',[tini:tfin],Mini,op);
C = M(:,1);
G = M(:,2);
E = M(:,3);
V = M(:,4);

Xc=(100*G).*V./(1.111*(Cini*Vini + Cin*Fin.*t)); %conversão


tempo=t/60; %(h)
%Figuras...
subplot (2,2,1)
plot(tempo,M(:,1),'-b')
grid on
xlabel('t(h)')
ylabel('Celulose (g/L)')
subplot (2,2,2)
plot(tempo,M(:,2),'-b')
grid on
xlabel('t(h)')
ylabel('Glicose (g/L)')
subplot (2,2,3)
plot(tempo,M(:,3),'-b')
grid on
xlabel('t(h)')
ylabel('Enzima (g/L)')
subplot (2,2,4)
plot(tempo,M(:,4),'-b')
grid on
xlabel('t(h)')
ylabel('Volume (L)')

plot (tempo,Xc,'-b')
grid on
xlabel('t(h)')
ylabel('Conversão (%)')


