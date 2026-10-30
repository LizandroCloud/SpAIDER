%Script to run degasser Level-1 MATLAB S-Function
clear

%Physical properties
global PP
PP=blackoil;

%Setting model parameters
par=cell(6,1);
par{1}=141.5/(30+131.5);                    %d60, Oil specific gravity []
par{2}=(1/1.8)*(11.8*par{1})^3;             %Tbn, Oil normal boiling point [K]
par{3}=255;     				  			%Cws, Water salt concentration [g/L]
par{4}=1.05;                                %dg, Gas specific gravity []
par{5}=8.00;                                %C, Vessel length [m]
par{6}=2.50;                                %D, Vessel diameter [m]

%Model states linearization point
hL=1.50;                 %Liquid level                       [m]
BSW=0.15;                %Basic sediments & water            []
P=801.325;               %Pressure                           [kPa]
T=80+273.15;             %Temperature                        [K]

%Model inputs linearization point
Qf_in=22000;              %Inlet liquid std. vol. flow rate  [Sm^3/d]
BSW_in=BSW;               %Inlet basic sediments & water     []
Qg_in=5*Qf_in*(1-BSW_in); %Inlet gas std. vol. flow rate     [Sm^3/d]
T_in=T;                   %Inlet temperature                 [K]
Qf=Qf_in;                 %Liquid std. vol. flow rate        [Sm^3/d]
Qg=Qg_in;                 %Gas std. vol. flow rate           [Sm^3/d]

%Getting state-space linear model
cinit=[hL; BSW; P; T];
u0=[Qf_in; Qg_in; BSW_in; T_in; Qf; Qg];
sys=sfunc2ss(@degasser,cinit,u0,par)