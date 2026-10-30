%Script to run uotr Level-1 MATLAB S-Function
clear

%Physical properties
global PP
PP=blackoil;

%----- Setting model parameters:
%black-oil physical properties
d60	= 141.5/(30+131.5);     		%Oil specific gravity 				[]
Tbn	= (1/1.8)*(11.8*d60)^3;			%Oil normal boiling point 			[K]
Cws	= [255,25];     				%Water salt concentration 			[g/L]
dg	= [0.7,1.05,1.5];       		%Gas specific gravity 				[]
	
%triphase_separator SG_001
SG_001.LCS		= 16.5;     %Separation chamber (CS) length 	[m]
SG_001.LCL		= 2.50;     %Lights chamber (CL) length 		[m]
SG_001.D		= 4.30;     %Vessel diameter 					[m]
SG_001.hweir	= 2.15;     %Weir height 						[m]
	
%valv_liq LV_001	
LV_001.Cv		= 1200;     %Valve constant US units 			[gpm]
LV_001.FL		= 0.91;     %Liquid pressure rec. factor		[]
	
%valv_water LV_002	
LV_002.Cv		= 500;      %Valve constant US units 			[gpm]
LV_002.FL		= 0.91;     %Liquid pressure rec. factor		[]
	
%valv_gas PV_001	
PV_001.Cv		= 1200;     %Valve constant US units 			[gpm]
PV_001.Fgm		= 0.80;     %Spec. heat cap. ratio factor 		[]
PV_001.DPrT		= 0.70;     %Pressure drop ratio factor 		[]
	
%hexch_liq P_001	
P_001.UA		= 1600;     %Global heat transfer coeff. 		[kW/K]
	
%degasser V_TO_001	
V_TO_001.C		= 9.45;     %Vessel length 						[m]
V_TO_001.D		= 3.15;     %Vessel diameter 					[m]

%oil_treater TO_001
TO_001.C		= 15.0;     %Vessel length                  	[m]
TO_001.D		= 4.30;     %Vessel diameter                    [m]
TO_001.tRL50	= 210;     	%Res. time for water half-removal	[s]
TO_001.gamma	= 3.6;     	%Water removal eff. correl. 2nd par []

%valv_liq LV_003	
LV_003.Cv		= 1800;     %Valve constant US units 			[gpm]
LV_003.FL		= 0.91;     %Liquid pressure rec. factor		[]

%valv_water LV_004	
LV_004.Cv		= 270;      %Valve constant US units 			[gpm]
LV_004.FL		= 0.91;     %Liquid pressure rec. factor		[]

%pump_water B_001	
B_001.h0		= 300;      %Shut-off head                      [m]
B_001.A			= 5.184e3;  %Pressure-flow pump relation coeff.	[s^2/m^5]

%valv_gas PV_002
PV_002.Cv		= 60;     	%Valve constant US units 			[gpm]
PV_002.Fgm		= 0.80;     %Spec. heat cap. ratio factor 		[]
PV_002.DPrT		= 0.70;     %Pressure drop ratio factor 		[]

%degasser V_TO_002	
V_TO_002.C		= 9.45;     %Vessel length 						[m]
V_TO_002.D		= 3.15;     %Vessel diameter 					[m]

%oil_treater TO_002
TO_002.C		= 15.0;     %Vessel length                  	[m]
TO_002.D		= 4.30;     %Vessel diameter                    [m]
TO_002.tRL50	= 210;     	%Res. time for water half-removal	[s]
TO_002.gamma	= 3.6;     	%Water removal eff. correl. 2nd par []

%valv_liq LV_005
LV_005.Cv		= 1800;     %Valve constant US units 			[gpm]
LV_005.FL		= 0.91;     %Liquid pressure rec. factor		[]

%valv_water LV_006	
LV_006.Cv		= 60;      %Valve constant US units 			[gpm]
LV_006.FL		= 0.91;     %Liquid pressure rec. factor		[]

%pump_water B_002	
B_002.h0		= 300;      %Shut-off head                      [m]
B_002.A			= 1.5552e5; %Pressure-flow pump relation coeff.	[s^2/m^5]

%pump_liq B_003	
B_003.h0		= 50;      	%Shut-off head                      [m]
B_003.A			= 1.5552e1;	%Pressure-flow pump relation coeff.	[s^2/m^5]
	
%valv_gas PV_003
PV_003.Cv		= 50;     	%Valve constant US units 			[gpm]
PV_003.Fgm		= 0.80;     %Spec. heat cap. ratio factor 		[]
PV_003.DPrT		= 0.70;     %Pressure drop ratio factor 		[]

%vapor recovery unit UVRC
UVRC.Nr_Max 		= 12000;	%Maximum shitf assemble rotation 	[rpm]
UVRC.nS 			= 2;		%Number of compression stages 		[]

%compressor UVRC.comp(k)
UVRC.comp(1).Fgm     	= 1.2;     			%Specific heat capacity ratio factor 	[]
UVRC.comp(1).PrCC		= [0.35,4.5,1.5];  	%Pres. ratio charac. curve par.   		[]
UVRC.comp(1).EffSPOL	= [0.35,1.1,7.0];   %Isentropic eff. polynomial par.   		[]

UVRC.comp(2).Fgm     	= 1.2;     			%Specific heat capacity ratio factor 	[]
UVRC.comp(2).PrCC		= [0.35,4.5,1.5];  	%Pres. ratio charac. curve par.   		[]
UVRC.comp(2).EffSPOL	= [0.35,1.1,7.0];   %Isentropic eff. polynomial par.   		[]

%valv_gas UVRC.ascv(k)
UVRC.ascv(1).Cv			= 300;              %Valve constant US units 				[gpm]
UVRC.ascv(1).Fgm		= 1.2;              %Spec. heat cap. ratio factor 			[]
UVRC.ascv(1).DPrT		= 0.70;     		%Pressure drop ratio factor 			[]
	
UVRC.ascv(2).Cv			= 300;              %Valve constant US units 				[gpm]
UVRC.ascv(2).Fgm		= 1.2;              %Spec. heat cap. ratio factor 			[]
UVRC.ascv(2).DPrT		= 0.70;     		%Pressure drop ratio factor 			[]

%Parameter structure
par.d60 		= d60;
par.Tbn 		= Tbn;
par.Cws 		= Cws;
par.dg 			= dg;
par.SG_001		= {SG_001.LCS, SG_001.LCL, SG_001.D, SG_001.hweir};
par.LV_001		= {LV_001.Cv, LV_001.FL};
par.LV_002		= {LV_002.Cv, LV_002.FL};
par.PV_001		= {PV_001.Cv, PV_001.Fgm, PV_001.DPrT};
par.P_001		= P_001.UA;
par.V_TO_001	= {V_TO_001.C, V_TO_001.D};
par.TO_001		= {TO_001.C, TO_001.D, TO_001.tRL50, TO_001.gamma};
par.LV_003		= {LV_003.Cv, LV_003.FL};
par.LV_004		= {LV_004.Cv, LV_004.FL};
par.B_001		= {B_001.h0, B_001.A};
par.PV_002		= {PV_002.Cv, PV_002.Fgm, PV_002.DPrT};
par.V_TO_002	= {V_TO_002.C, V_TO_002.D};
par.TO_002		= {TO_002.C, TO_002.D, TO_002.tRL50, TO_002.gamma};
par.LV_005		= {LV_005.Cv, LV_005.FL};
par.LV_006		= {LV_006.Cv, LV_006.FL};
par.B_002		= {B_002.h0, B_002.A};
par.B_003		= {B_003.h0, B_003.A};
par.PV_003		= {PV_003.Cv, PV_003.Fgm, PV_003.DPrT};
par.UVRC.Nr_Max	= UVRC.Nr_Max;
par.UVRC.nS		= UVRC.nS;
par.UVRC.comp   = UVRC.comp;
par.UVRC.ascv   = UVRC.ascv;

%----- States linearization point:
%triphase_separator SG_001
SG_001.hWCS   	= 0.15; 	%Water phase level in CS        	[m]
SG_001.BSWLCS 	= 0.20;		%Basic sediments & water in CS 		[]
SG_001.hCS = 1.001*SG_001.hweir;%Liquid level in CS 			[]
SG_001.hCL    	= 1.5;		%Liquid level in CL             	[m]
SG_001.BSWLCL	= 0.20;		%Basic sediments & water in CL 		[]
SG_001.P       	= 2101.325;	%Pressure                       	[kPa]
SG_001.T       	= 300;		%Temperature                    	[K]
SG_001.x0 = [...	
	SG_001.hWCS; ...
	SG_001.BSWLCS; ... 
	SG_001.hCS; ...
	SG_001.hCL; ...
	SG_001.BSWLCL; ...
	SG_001.P; ...
	SG_001.T];
	
%degasser V_TO_001
V_TO_001.hL   	= 1.00; 	%Liquid level 				       	[m]
V_TO_001.BSW 	= 0.20;		%Basic sediments & water			[]
V_TO_001.P      = 801.325;	%Pressure                       	[kPa]
V_TO_001.T      = 360;		%Temperature                    	[K]
V_TO_001.x0 = [...
	V_TO_001.hL; ...
	V_TO_001.BSW; ...
	V_TO_001.P; ...
	V_TO_001.T];

%oil_treater TO_001
TO_001.hW   	= 2.00; 	%Water phase level 				    [m]
TO_001.T      	= 360;		%Temperature                    	[K]
TO_001.x0 = [...	
	TO_001.hW; ...
	TO_001.T];

%degasser V_TO_002
V_TO_002.hL   	= 1.00; 	%Liquid level 				       	[m]
V_TO_002.BSW 	= 0.05;		%Basic sediments & water			[]
V_TO_002.P      = 301.325;	%Pressure                       	[kPa]
V_TO_002.T      = 360;		%Temperature                    	[K]
V_TO_002.x0 = [...
	V_TO_002.hL; ...
	V_TO_002.BSW; ...
	V_TO_002.P; ...
	V_TO_002.T];

%oil_treater TO_002
TO_002.hW   	= 2.00; 	%Water phase level 				    [m]
TO_002.T      	= 360;		%Temperature                    	[K]
TO_002.x0 = [...	
	TO_002.hW; ...
	TO_002.T];

%----- Inputs linearization point:
%triphase_separator SG_001 feed stream
feed.Qf_in		= 24000;    %Inlet liquid std. vol. flow rate  	[Sm^3/d]
feed.BSW_in		= 0.1;      %Inlet basic sediments & water     	[]
feed.Qg_in 		= 250*feed.Qf_in*(1-feed.BSW_in);
							%Inlet gas std. vol. flow rate 		[Sm^3/d]
feed.T_in		= 290;		%Inlet temperature                 	[K]

%valv_liq LV_001
LV_001.c		= 0.5;		%Input signal from controller      	[]	

%valv_water LV_002
LV_002.P_out	= 801.325;  %Outlet pressure                   	[kPa]
LV_002.c		= 0.5;		%Input signal from controller      	[]

%valv_gas PV_001
%PV_001.P_out	= 2001.325;	%Outlet pressure                   	[kPa]
PV_001.P_out	= 501.325;	%Outlet pressure                   	[kPa]
PV_001.c		= 0.7;		%Input signal from controller      	[]

%hexch_liq P_001
P_001.DP_C  	= 120;		%Cold side pressure drop			[kPa]
P_001.DP_H  	= 120;		%Cold side pressure drop			[kPa]

%heater_liq P_002
P_002.DP  		= 120;		%Production heater pressure drop	[kPa]
P_002.Q_in    	= 15e3;		%Production heater added heat       [kW]

%valv_liq LV_003
LV_003.c     	= 0.5;		%Input signal from controller      	[]

%valv_water LV_004
LV_004.c     	= 0.5;		%Input signal from controller      	[]

%valv_gas PV_002
PV_002.P_out 	= 101.325;	%Outlet pressure                   	[kPa]
%PV_002.c     	= 0.5;		%Input signal from controller      	[]
PV_002.c     	= 1e-5;		%Input signal from controller      	[]

%mixer_W3ph D_001 desalt water stream
dslt.Qw 		= 1200;		%Std. vol. flow rate 				[Sm^3/d]
dslt.T 			= 360;		%Temperature 						[K]

%valv_liq LV_005
LV_005.c     	= 0.5;		%Input signal from controller      	[]

%valv_water LV_006
LV_006.c     	= 0.5;		%Input signal from controller      	[]

%Stock tank stream
stk.P = 101.325;			%Stock tank pressure 				[kPa]

%valv_gas PV_003
PV_003.P_out 	= 101.325;	%Outlet pressure                   	[kPa]
%PV_003.c     	= 0.5;		%Input signal from controller      	[]
PV_003.c     	= 1e-5;		%Input signal from controller      	[]

%vapor recovery unit UVRC
UVRC.MV_in 			= 0.75;		%Input signal from capacity control	[]

%compressor UVRC.cool(k)
UVRC.cool(1).Q_out 	= 300;		%Heat removed from system 			[kW]
UVRC.cool(1).DP 	= 30;		%Pressure drop 						[kPa]
UVRC.cool(2).Q_out 	= 650;		%Heat removed from system 			[kW]
UVRC.cool(2).DP 	= 30;		%Pressure drop 						[kPa]

%valv_gas UVRC.ascv(k)
UVRC.ascv(1).c 		= 0.1;		%Valve stem position 				[]
UVRC.ascv(2).c 		= 0.1;		%Valve stem position 				[]

%cooler_gas P_003
P_003.DP			= 30;		%Pressure drop 						[kPa]
P_003.Q_out			= 1200;		%Heat removed from system 			[kW]

%----- Outputs linearization point:
%heater_liq P_002
P_002.T  		= 360;		%Temperature 						[K]

%cooler_gas P_003
P_003.T  		= 310;		%Temperature 						[K]

%compressor UVRC.cool(k)
UVRC.cool(1).T 	= 320;		%Temperature 						[K]
UVRC.cool(2).T 	= 320;		%Temperature 						[K]

%----- Steady-state computation:
x0 = [...
	SG_001.x0; ...
	V_TO_001.x0; ...
	TO_001.x0; ...
	V_TO_002.x0; ...
	TO_002.x0];

u0 = [...
	feed.Qf_in; feed.Qg_in; feed.BSW_in; feed.T_in; ...
	LV_001.c; ...
	LV_002.P_out; LV_002.c; ...
	PV_001.P_out; PV_001.c; ...
	stk.P; ...
	P_001.DP_C; P_001.DP_H; ...
	P_002.DP; P_002.Q_in; ...
	LV_003.c; ...
	LV_004.c; ...
	PV_002.P_out; PV_002.c; ...
	dslt.Qw; dslt.T; ...
	LV_005.c; ...
	LV_006.c; ...
	PV_003.P_out; PV_003.c; ...
	UVRC.cool(1).Q_out; UVRC.cool(1).DP; ...
	UVRC.cool(2).Q_out; UVRC.cool(2).DP; ...
	UVRC.ascv(1).c; ...
	UVRC.ascv(2).c; ...
	UVRC.MV_in; ...
	P_003.Q_out; P_003.DP];

yss = [...
	P_002.T; ...
	P_003.T; ...
	UVRC.cool(1).T; ...
	UVRC.cool(2).T];

[x0,u0]=uotr_ss(x0,u0,yss,par);

%----- State-space linear model:
sys=sfunc2ss(@uotr,x0,u0,par)