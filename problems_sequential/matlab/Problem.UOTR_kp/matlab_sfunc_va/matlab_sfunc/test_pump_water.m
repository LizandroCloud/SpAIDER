%Script to run valv_liq Level-1 MATLAB S-Function
clear

%Physical properties
global PP
PP=blackoil;

%Setting model parameters
par=cell(3,1);
par{1}=255;     				  %Cws, Water salt concentration [g/L]
par{2}=220;                       %h0, pump head [m]
par{3}=9e-4*3600^2;               %A, Pressure-flow pump relation coefficient [s^2/m^5]

%Model inputs linearization point
Qw_in=5000;			%Inlet water std. vol. flow rate   [Sm^3/d]
P_in=801.325;       %Inlet pressure                    [kPa]
T_in=360;           %Inlet temperature                 [K]

%Getting state-space linear model
cinit=[];
u0=[Qw_in; P_in; T_in];
sys=sfunc2ss(@pump_water,cinit,u0,par)