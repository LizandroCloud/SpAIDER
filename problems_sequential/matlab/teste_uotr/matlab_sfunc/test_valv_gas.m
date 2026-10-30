%Script to run valv_gas Level-1 MATLAB S-Function
clear

%Physical properties
global PP
PP=blackoil;

%Setting model parameters
par=cell(5,1);
par{1}=0.7;                       %dg, Gas specific gravity []
par{2}=1500;                      %Cv, Valve constant US units [gpm]
par{3}=0.80;                      %Fgm, Specific heat capacity ratio factor []
par{4}=0.70;                      %DPrT, Pressure drop ratio factor []

%Model inputs linearization point
P_in=2101.325;      		%Inlet pressure                    [kPa]
T_in=300;           		%Inlet temperature                 [K]
P_out=2000.325;     		%Outlet pressure                   [kPa]
c=0.5;					%Input signal from controller      []

%Getting state-space linear model
u0=[P_in; T_in; P_out; c];
sys=sfunc2ss(@valv_gas,[],u0,par)