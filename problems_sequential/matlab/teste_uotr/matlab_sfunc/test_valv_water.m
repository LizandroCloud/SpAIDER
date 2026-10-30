%Script to run valv_water Level-1 MATLAB S-Function
clear

%Physical properties
global PP
PP=blackoil;

%Setting model parameters
par=cell(4,1);
par{1}=255;     				  %Cws, Water salt concentration [g/L]
par{2}=500;                       %Cv, Valve constant US units [gpm]
par{3}=0.91;                      %FL, Liquid pressure recovery factor     []

%Model inputs linearization point
P_in=2101.325;      	%Inlet pressure                    [kPa]
T_in=300;           	%Inlet temperature                 [K]
P_out=801.325;      	%Outlet pressure                   [kPa]
c=0.5;					%Input signal from controller      []

%Getting state-space linear model
u0=[P_in; T_in; P_out; c];
sys=sfunc2ss(@valv_water,[],u0,par)