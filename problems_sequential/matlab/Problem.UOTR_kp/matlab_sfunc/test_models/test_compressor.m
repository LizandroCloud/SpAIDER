%Script to run compressor Level-1 MATLAB S-Function
clear

%Physical properties
global PP
PP=blackoil;

%Setting model parameters
par=cell(5,1);
par{1}=0.7;                 %dg, Gas specific gravity []
par{2}=0.80;                %Fgm, Specific heat capacity ratio factor []
par{3}=12000;               %Nr_Max, Maximum compressor rotation [rpm]
par{4}=[0.35,1.3,0.8];      %PrCC, Pres. ratio charac. curve par. []
par{5}=[0.35,0.36,0.6];     %EffSPOL, Isentropic eff. polynomial par. []

%Model inputs linearization point
P_in=2101.325;      		%Inlet pressure                    [kPa]
T_in=300;           		%Inlet temperature                 [K]
P_out=8101.325;     		%Outlet pressure                   [kPa]
MV_in=0.9;					%Input signal from controller      []

%Getting state-space linear model
u0=[P_in; T_in; P_out; MV_in];
sys=sfunc2ss(@compressor,[],u0,par)