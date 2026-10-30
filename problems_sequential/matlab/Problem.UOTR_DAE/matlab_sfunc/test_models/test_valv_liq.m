%Script to run valv_liq Level-1 MATLAB S-Function
clear

%Physical properties
global PP
PP=blackoil;

%Setting model parameters
par=cell(7,1);
par{1}=141.5/(30+131.5);          %d60, Oil specific gravity []
par{2}=(1/1.8)*(11.8*par{1})^3;   %Tbn, Oil normal boiling point [K]
par{3}=255;     				  %Cws, Water salt concentration [g/L]
par{4}=0.7;                       %dg, Gas specific gravity []
par{5}=900;                       %Cv, Valve constant US units [gpm]
par{6}=0.91;                      %FL, Liquid pressure recovery factor     []

%Model inputs linearization point
BSW_in=0.20;        %Inlet basic sediments & water     []
P_in=2101.325;      %Inlet pressure                    [kPa]
T_in=300;           %Inlet temperature                 [K]
P_out=801.325;      %Outlet pressure                   [kPa]
c=0.5;				%Input signal from controller      []

%Getting state-space linear model
u0=[BSW_in; P_in; T_in; P_out; c];
sys=sfunc2ss(@valv_liq,[],u0,par)