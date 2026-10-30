% MODELO PSEUDO-HOMOGÊNEO HIDRÓLISE AVICEL
% Modo:Batelada alimentada

function dM = Avicel(t,M)

% DADOS ENTRADA
Fenz = 0.0001;      % Vazão de enzima (L/min)
Fsol = 0.0004;      % Vazão de sólidos (L/min)
Fout = 0.0;         % Vazão de saida do reator (L/min)
Eenz = 1.937*(5/1); % [] enzima entra mixer (g/L)
Csol = 80*(5/4);    % [] celulose entra mixer (g/L)
Gsol = 0*(5/4);     % [] glicose entra mixer (g/L)

% PARÂMETROS
k = 0.0257458;  %(min^-1)
Km = 11.6035;   %(g/L)
Kig = 5.38839;  %(g/L)

% VARIÁVEIS
C = M(1); %[] Substrato-celulose (g/L)
G = M(2); %[] Glicose (g/L)
E = M(3); %[] Enzima (g/L)
V = M(4); %Volume (L)

% EQUAÇÕES
% Eq. Mixer:
Fin = Fsol + Fenz;
Cin = (Fsol/Fin)*Csol;
Gin = (Fsol/Fin)*Gsol;
Ein = (Fenz/Fin)*Eenz;
% Eq. Reator:
% C->G. (inibição competitiva da G)
r= k*E*C/(Km*(1+ G/Kig) + C);       %(g/L*min)  
%EDO's
dM(1)= (Fin/V)*(Cin-C) - r; %C
%dM(1)=  - r; %C
dM(2)= (Fin/V)*(Gin-G) + 1.111*r; %G
%dM(2)= + 1.111*r; %G
dM(3)= (Fin/V)*(Ein-E);     %E
dM(4)= (Fin-Fout);          %V
dM   = dM';

end