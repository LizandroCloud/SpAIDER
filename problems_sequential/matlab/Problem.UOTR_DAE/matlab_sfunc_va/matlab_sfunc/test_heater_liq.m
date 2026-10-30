%Script to run heater_liq Level-1 MATLAB S-Function
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

%Model inputs linearization point
Qf_in=24000;              	%Inlet liquid std. vol. flow rate  [Sm^3/d]
BSW_in=0.2;               	%Inlet basic sediments & water     []
Qg_in=5*Qf_in*(1-BSW_in); 	%Inlet gas std. vol. flow rate     [Sm^3/d]
P_in=2101.325;      		%Inlet pressure                    [kPa]
T_in=300;           		%Inlet temperature                 [K]
P_out=801.325;      		%Outlet pressure                   [kPa]
Q_in=20e3					%Heat applied to the system        [kW]

%Getting state-space linear model
u0=[Qf_in; Qg_in; BSW_in; P_in; T_in; P_out; Q_in];
sys=sfunc2ss(@heater_liq,[],u0,par)