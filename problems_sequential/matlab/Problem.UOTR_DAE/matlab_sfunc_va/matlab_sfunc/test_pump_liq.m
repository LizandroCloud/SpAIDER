%Script to run valv_liq Level-1 MATLAB S-Function
clear

%Physical properties
global PP
PP=blackoil;

%Setting model parameters
par=cell(6,1);
par{1}=141.5/(30+131.5);          %d60, Oil specific gravity []
par{2}=(1/1.8)*(11.8*par{1})^3;   %Tbn, Oil normal boiling point [K]
par{3}=255;     				  %Cws, Water salt concentration [g/L]
par{4}=0.7;                       %dg, Gas specific gravity []
par{5}=220;                       %h0, pump head [m]
par{6}=9e-4*3600^2;               %A, Pressure-flow pump relation coefficient [s^2/m^5]

%Model inputs linearization point
Qf_in=22000;		%Inlet liq. std. vol. flow rate 	[Sm^3/d]
BSW_in=0.25;        %Inlet basic sediments & water     []
P_in=801.325;       %Inlet pressure                    [kPa]
T_in=300;           %Inlet temperature                 [K]

%Getting state-space linear model
cinit=[];
u0=[Qf_in; BSW_in; P_in; T_in];
sys=sfunc2ss(@pump_liq,cinit,u0,par)