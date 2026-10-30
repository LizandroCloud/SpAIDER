%Script to run oil_treater Level-1 MATLAB S-Function
clear

%Physical properties
global PP
PP=blackoil;

%Setting model parameters
par=cell(8,1);
par{1}=141.5/(30+131.5);          %d60, Oil specific gravity []
par{2}=(1/1.8)*(11.8*par{1})^3;   %Tbn, Oil normal boiling point [K]
par{3}=255;     				  %Cws, Water salt concentration [g/L]
par{4}=1.05;                      %dg, Gas specific gravity []
par{5}=15.0;                      %C, vessel length [m]
par{6}=4.00;                      %D, Vessel diameter [m]
par{7}=210;                       %tRL50, Residence time for water half-removal  [s]
par{8}=3.6;                       %gamma, Water removal eff. correl. 2nd parameter []

%Model states linearization point
hW   	=   2.0; 			%Water phase level        			[m]
T       =   360;			%Temperature                        [K]

%Model inputs linearization point
BSW_in=0.2;               	%Inlet basic sediments & water     	[]
P_in=801.325;				%Inlet pressure 				 	[kPa]
T_in=T;                   	%Inlet temperature                 	[K]
Qf=22000;              		%Out. liquid std. vol. flow rate   	[Sm^3/d]
BSW=0.005;
effWL=1-BSW/BSW_in;
Qw=Qf*BSW_in*effWL;          %Out. water std. vol. flow rate   	[Sm^3/d]

%Getting state-space linear model
cinit=[hW; T];
u0=[BSW_in; P_in; T_in; Qf; Qw];
sys=sfunc2ss(@oil_treater,cinit,u0,par)