%Script to run hexch_liq Level-1 MATLAB S-Function
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
par{5}=1600;                      %UA, Global heat transfer coeff. [kW/K]

%Model inputs linearization point
Qf_inC=20000;              	%Cold Inlet liquid std. vol. flow rate  [Sm^3/d]
BSW_inC=0.2;               	%Cold Inlet basic sediments & water     []
Rgo_inC=5;
Qg_inC=Rgo_inC*Qf_inC*(1-BSW_inC); 	%Cold Inlet gas std. vol. flow rate     [Sm^3/d]
P_inC=801.325;      		%Cold Inlet pressure                    [kPa]
T_inC=310;           		%Cold Inlet temperature                 [K]
P_outC=P_inC-120;      		%Cold Outlet pressure                   [kPa]

Qf_inH=16000;              	%Hot Inlet liquid std. vol. flow rate  [Sm^3/d]
BSW_inH=0.005;              %Hot Inlet basic sediments & water     []
Rgo_inH=1e-6;
Qg_inH=Rgo_inH*Qf_inH*(1-BSW_inH); 	%Hot Inlet gas std. vol. flow rate     [Sm^3/d]
P_inH=301.325;      		%Hot Inlet pressure                    [kPa]
T_inH=360;           		%Hot Inlet temperature                 [K]
P_outH=P_inH-120;      		%Hot Outlet pressure                   [kPa]

%Getting state-space linear model
u0=[Qf_inC; Qg_inC; BSW_inC; P_inC; T_inC; P_outC; Qf_inH; Qg_inH; BSW_inH; P_inH; T_inH; P_outH];
sys=sfunc2ss(@hexch_liq,[],u0,par)