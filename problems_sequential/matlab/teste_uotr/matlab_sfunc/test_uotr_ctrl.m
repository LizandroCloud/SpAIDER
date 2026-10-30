%Script to run uotr_ctrl Level-1 MATLAB S-Function
clear, close all

%Physical properties
global PP
PP=blackoil;
%{
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
%{
%valv_gas PV_002
PV_002.Cv		= 60;     	%Valve constant US units 			[gpm]
PV_002.Fgm		= 0.80;     %Spec. heat cap. ratio factor 		[]
PV_002.DPrT		= 0.70;     %Pressure drop ratio factor 		[]
%}
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
%{
%valv_gas PV_003
PV_003.Cv		= 50;     	%Valve constant US units 			[gpm]
PV_003.Fgm		= 0.80;     %Spec. heat cap. ratio factor 		[]
PV_003.DPrT		= 0.70;     %Pressure drop ratio factor 		[]
%}
%vapor recovery unit UVRC
UVRC.Nr_Max 	= 12000;		%Maximum shitf assemble rotation 	[rpm]
UVRC.nS 		= 2;			%Number of compression stages 		[]
UVRC.Qout_Max	= [2e3, 2e3];	%Max heat removed from system 		[kW]

%scrubber_gas UVRC.scrb(k)
UVRC.scrb(1).C 		= 3.0;				%Scrubber length 						[m]
UVRC.scrb(1).D 		= 1.25;				%Scrubber diameter 						[m]

UVRC.scrb(2).C 		= 3.0;				%Scrubber length 						[m]
UVRC.scrb(2).D 		= 1.25;				%Scrubber diameter 						[m]

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

%scrubber_gas V_001
V_001.C 				= 6.00;				%Vessel length 							[m]
V_001.D 				= 2.50;				%Vessel diameter 						[m]

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
%par.PV_002		= {PV_002.Cv, PV_002.Fgm, PV_002.DPrT};
par.V_TO_002	= {V_TO_002.C, V_TO_002.D};
par.TO_002		= {TO_002.C, TO_002.D, TO_002.tRL50, TO_002.gamma};
par.LV_005		= {LV_005.Cv, LV_005.FL};
par.LV_006		= {LV_006.Cv, LV_006.FL};
par.B_002		= {B_002.h0, B_002.A};
par.B_003		= {B_003.h0, B_003.A};
%par.PV_003		= {PV_003.Cv, PV_003.Fgm, PV_003.DPrT};
par.UVRC.Nr_Max	= UVRC.Nr_Max;
par.UVRC.Qout_Max = UVRC.Qout_Max;
par.UVRC.nS		= UVRC.nS;
par.UVRC.scrb   = UVRC.scrb;
par.UVRC.comp   = UVRC.comp;
par.UVRC.ascv   = UVRC.ascv;
par.V_001 		= {V_001.C, V_001.D};
%{
%----- Linearization points:
%scrubber_gas V_001
V_001.P 		= 2101.325;				%Pressure 		[kPa]
V_001.T 		= 310;          		%Temperature    [K]
V_001.x0 = [...
	V_001.P; ...
	V_001.T];

%cooler_gas P_003
P_003.DP			= 30;		%Pressure drop 						[kPa]
P_003.Q_out			= 1.2e3;	%Heat removed from system 			[kW]

%PID_pos TIC_002
TIC_002.MVinput			= 1.2/2.5;
TIC_002.derivFilt 		= eps;
TIC_002.AWBTintterm 	= -TIC_002.MVinput; %TIC_002.AWBTintterm 	= eps;
TIC_002.x0 = [TIC_002.AWBTintterm; TIC_002.derivFilt; TIC_002.MVinput];

%triphase_separator SG_001
SG_001.hWCS   	= 0.15; 				%Water phase level in CS        	[m]
SG_001.BSWLCS 	= 0.20;					%Basic sediments & water in CS 		[]
SG_001.hCS = 1.001*SG_001.hweir;		%Liquid level in CS 			[]
SG_001.hCL    	= 1.5;					%Liquid level in CL             	[m]
SG_001.BSWLCL	= 0.20;					%Basic sediments & water in CL 		[]
SG_001.P 	= V_001.P + P_003.DP + 1e-3;%Pressure                       	[kPa]
SG_001.T       	= 310;					%Temperature                    	[K]
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

%compressor UVRC.cool(k)
UVRC.cool(1).Q_out 	= 0.4e3;	%Heat removed from system 			[kW]
UVRC.cool(1).DP 	= 30;		%Pressure drop 						[kPa]
UVRC.cool(2).Q_out 	= 1.25e3;	%Heat removed from system 			[kW]
UVRC.cool(2).DP 	= 30;		%Pressure drop 						[kPa]

%vapor recovery unit UVRC
UVRC.scrb(1).P 	= V_TO_002.P - UVRC.cool(1).DP - 1e-4;	%Pressure       [kPa]
UVRC.scrb(1).T 	= 320;                          %Temperature    [K]
UVRC.scrb(2).P  = V_TO_001.P - UVRC.cool(2).DP - 1e-4; %Pressure       [kPa]
UVRC.scrb(2).T  = 320;                          %Temperature    [K]
	
%P_pos UVRC.asc(k)
UVRC.asc(1).MVinput 	= 0.2;
UVRC.asc(2).MVinput 	= 0.2;

%PID_pos UVRC.TIC(k)
UVRC.TIC(1).MVinput		= UVRC.cool(1).Q_out/UVRC.Qout_Max(1);
UVRC.TIC(1).derivFilt 	= eps;
UVRC.TIC(1).AWBTintterm = -UVRC.TIC(1).MVinput; %UVRC.TIC(1).AWBTintterm = eps;

UVRC.TIC(2).MVinput		= UVRC.cool(2).Q_out/UVRC.Qout_Max(2);
UVRC.TIC(2).derivFilt 	= eps;
UVRC.TIC(2).AWBTintterm = -UVRC.TIC(2).MVinput; %UVRC.TIC(2).AWBTintterm = eps;

%PID_pos UVRC.MCC
UVRC.MCC.MVinput		= 0.75;
UVRC.MCC.derivFilt 		= eps;
UVRC.MCC.AWBTintterm 	= -UVRC.MCC.MVinput; %UVRC.MCC.AWBTintterm 	= eps;

UVRC.x0 = [...
	UVRC.scrb(1).P; UVRC.scrb(1).T; ...
	UVRC.scrb(2).P; UVRC.scrb(2).T; ...
	UVRC.asc(1).MVinput; UVRC.asc(2).MVinput; ...
	UVRC.TIC(1).AWBTintterm; UVRC.TIC(1).derivFilt; UVRC.TIC(1).MVinput; ...
	UVRC.TIC(2).AWBTintterm; UVRC.TIC(2).derivFilt; UVRC.TIC(2).MVinput; ...
	UVRC.MCC.AWBTintterm; UVRC.MCC.derivFilt; UVRC.MCC.MVinput ...
	];
	
%triphase_separator SG_001 feed stream
feed.Qf_in		= 24000;    %Inlet liquid std. vol. flow rate  	[Sm^3/d]
feed.BSW_in		= 0.1;      %Inlet basic sediments & water     	[]
feed.Rgo 		= 250;		%Gas-oil ratio 						[Sm^3/Sm^3]
							%Inlet gas std. vol. flow rate 		[Sm^3/d]
feed.T_in		= 290;		%Inlet temperature                 	[K]

%valv_liq LV_001
LV_001.c		= 0.5;		%Input signal from controller      	[]

%PID_pos LIC_001
LIC_001.MVinput			= LV_001.c;
LIC_001.derivFilt 		= eps;
LIC_001.AWBTintterm 	= -LIC_001.MVinput; %LIC_001.AWBTintterm 	= eps;
LIC_001.x0 = [LIC_001.AWBTintterm; LIC_001.derivFilt; LIC_001.MVinput];

%valv_water LV_002
LV_002.c		= 0.5;		%Input signal from controller      	[]
LV_002.P_out	= 801.325;  %Outlet pressure                   	[kPa]

%PID_pos LIC_002
LIC_002.MVinput			= LV_002.c;
LIC_002.derivFilt 		= eps;
LIC_002.AWBTintterm 	= -LIC_002.MVinput; %LIC_002.AWBTintterm 	= eps;
LIC_002.x0 = [LIC_002.AWBTintterm; LIC_002.derivFilt; LIC_002.MVinput];

%valv_gas PV_001
PV_001.P_out	= 501.325;	%Outlet pressure                   	[kPa]
PV_001.c		= 0.5;		%Input signal from controller      	[]

%PID_pos PIC_001
PIC_001.MVinput			= PV_001.c;
PIC_001.derivFilt 		= eps;
PIC_001.AWBTintterm 	= -PIC_001.MVinput; %PIC_001.AWBTintterm 	= eps;
PIC_001.x0 = [PIC_001.AWBTintterm; PIC_001.derivFilt; PIC_001.MVinput];

%hexch_liq P_001
P_001.DP_C  	= 120;		%Cold side pressure drop			[kPa]
P_001.DP_H  	= 120;		%Cold side pressure drop			[kPa]

%heater_liq P_002
P_002.DP  		= 120;		%Production heater pressure drop	[kPa]
P_002.Q_in    	= 15e3;		%Production heater added heat       [kW]

%PID_pos TIC_001
TIC_001.MVinput			= (P_002.Q_in-5e3)/(25e3-5e3);
TIC_001.derivFilt 		= eps;
TIC_001.AWBTintterm 	= TIC_001.MVinput; %TIC_001.AWBTintterm 	= eps;
TIC_001.x0 = [TIC_001.AWBTintterm; TIC_001.derivFilt; TIC_001.MVinput];

%valv_liq LV_003
LV_003.c     	= 0.5;		%Input signal from controller      	[]

%PID_pos LIC_003
LIC_003.MVinput			= LV_003.c;
LIC_003.derivFilt 		= eps;
LIC_003.AWBTintterm 	= -LIC_003.MVinput; %LIC_003.AWBTintterm 	= eps;
LIC_003.x0 = [LIC_003.AWBTintterm; LIC_003.derivFilt; LIC_003.MVinput];

%valv_water LV_004
LV_004.c     	= 0.5;		%Input signal from controller      	[]

%PID_pos LIC_004
LIC_004.MVinput			= LV_004.c;
LIC_004.derivFilt 		= eps;
LIC_004.AWBTintterm 	= -LIC_004.MVinput; %LIC_004.AWBTintterm 	= eps;
LIC_004.x0 = [LIC_004.AWBTintterm; LIC_004.derivFilt; LIC_004.MVinput];

%mixer_W3ph D_001 desalt water stream
dslt.Qw 		= 1200;		%Std. vol. flow rate 				[Sm^3/d]
dslt.T 			= 360;		%Temperature 						[K]

%valv_liq LV_005
LV_005.c     	= 0.5;		%Input signal from controller      	[]

%PID_pos LIC_005
LIC_005.MVinput			= LV_005.c;
LIC_005.derivFilt 		= eps;
LIC_005.AWBTintterm 	= -LIC_005.MVinput; %LIC_005.AWBTintterm 	= eps;
LIC_005.x0 = [LIC_005.AWBTintterm; LIC_005.derivFilt; LIC_005.MVinput];

%valv_water LV_006
LV_006.c     	= 0.5;		%Input signal from controller      	[]

%PID_pos LIC_006
LIC_006.MVinput			= LV_006.c;
LIC_006.derivFilt 		= eps;
LIC_006.AWBTintterm 	= -LIC_006.MVinput; %LIC_006.AWBTintterm 	= eps;
LIC_006.x0 = [LIC_006.AWBTintterm; LIC_006.derivFilt; LIC_006.MVinput];

%Stock tank stream
stk.P = 101.325;			%Stock tank pressure 				[kPa]

%----- Steady-state computation:
x0 = [...
	SG_001.x0; ...
	V_TO_001.x0; ...
	TO_001.x0; ...
	V_TO_002.x0; ...
	TO_002.x0; ...
	UVRC.x0; ...
	V_001.x0; ...
	LIC_001.x0; ...
	LIC_002.x0; ...
	LIC_003.x0; ...
	LIC_004.x0; ...
	LIC_005.x0; ...
	LIC_006.x0; ...
	PIC_001.x0; ...
	TIC_001.x0; ...
	TIC_002.x0  ...
	];

u0 = [...
	feed.Qf_in; feed.Rgo; feed.BSW_in; feed.T_in; ...
	LV_002.P_out;
	PV_001.P_out;
	stk.P; ...
	P_001.DP_C; P_001.DP_H; ...
	P_002.DP; ...
	dslt.Qw; dslt.T; ...
	UVRC.cool(1).DP; ...
	UVRC.cool(2).DP; ...
	P_003.DP ...
	];

yss=[];

[x0,u0]=uotr_ctrl_ss(x0,u0,yss,par);
x0=[x0;eps];	%Adding state Vg,flare as null in t0
%}

%----- Linearization points (from uotr_ss.mat):
load('uotr_ss.mat','x0','u0');

%States:
%triphase_separator SG_001
SG_001.x0 			= x0(1:7);
SG_001.hWCS  		= SG_001.x0(1);
SG_001.BSWLCS 		= SG_001.x0(2);
SG_001.hCS  		= SG_001.x0(3);
SG_001.hCL  		= SG_001.x0(4);
SG_001.BSWLCL 		= SG_001.x0(5);
SG_001.P       		= SG_001.x0(6);
SG_001.T			= SG_001.x0(7);

%degasser V_TO_001
V_TO_001.x0			= x0(8:11);
V_TO_001.hL 		= V_TO_001.x0(1);
V_TO_001.BSW		= V_TO_001.x0(2);
V_TO_001.P  		= V_TO_001.x0(3);
V_TO_001.T  		= V_TO_001.x0(4);

%oil_treater TO_001
TO_001.x0			= x0(12:13);
TO_001.hW 			= TO_001.x0(1);
TO_001.T  			= TO_001.x0(2);

%degasser V_TO_002
V_TO_002.x0			= x0(14:17);
V_TO_002.hL 		= V_TO_002.x0(1);
V_TO_002.BSW		= V_TO_002.x0(2);
V_TO_002.P  		= V_TO_002.x0(3);
V_TO_002.T  		= V_TO_002.x0(4);

%oil_treater TO_002
TO_002.x0			= x0(18:19);
TO_002.hW 			= TO_002.x0(1);
TO_002.T  			= TO_002.x0(2);

%vapor recovery unit UVRC
UVRC.x0					= x0(20:25);
UVRC.scrb(1).P 			= UVRC.x0(1);
UVRC.scrb(1).T 			= UVRC.x0(2);
UVRC.scrb(2).P 			= UVRC.x0(3);
UVRC.scrb(2).T 			= UVRC.x0(4);
UVRC.asc(1).MVinput 	= UVRC.x0(5);
UVRC.asc(2).MVinput 	= UVRC.x0(6);

%scrubber_gas V_001
V_001.x0			= x0(26:27);
V_001.P 			= V_001.x0(1);
V_001.T 			= V_001.x0(2);

%Inputs:
feed.Qf_in			= u0(1);
feed.Rgo			= u0(2);
feed.BSW_in			= u0(3);
feed.T_in			= u0(4);
LV_001.c			= u0(5);
LV_002.P_out		= u0(6);
LV_002.c			= u0(7);
PV_001.P_out		= u0(8);
PV_001.c			= u0(9);
stk.P  				= u0(10);
P_001.DP_C  		= u0(11);
P_001.DP_H  		= u0(12);
P_002.DP  			= u0(13);
P_002.Q_in    		= u0(14);
LV_003.c     		= u0(15);
LV_004.c     		= u0(16);
dslt.Qw    			= u0(17);
dslt.T     			= u0(18);
LV_005.c        	= u0(19);
LV_006.c        	= u0(20);
UVRC.cool(1).Q_out 	= u0(21);
UVRC.cool(1).DP   	= u0(22);
UVRC.cool(2).Q_out  = u0(23);
UVRC.cool(2).DP   	= u0(24);
UVRC.MV_in     		= u0(25);
P_003.Q_out			= u0(26);
P_003.DP			= u0(27);

%Control loops
%P_pos UVRC.asc(k)
%Already set in states section vapor recovery unit UVRC

%PI_pos UVRC.TIC(k)
UVRC.TIC(1).MVinput		= UVRC.cool(1).Q_out/UVRC.Qout_Max(1);
UVRC.TIC(1).derivFilt 	= eps;
UVRC.TIC(1).AWBTintterm = -UVRC.TIC(1).MVinput; %UVRC.TIC(1).AWBTintterm = eps;

UVRC.TIC(2).MVinput		= UVRC.cool(2).Q_out/UVRC.Qout_Max(2);
UVRC.TIC(2).derivFilt 	= eps;
UVRC.TIC(2).AWBTintterm = -UVRC.TIC(2).MVinput; %UVRC.TIC(2).AWBTintterm = eps;

%PI_pos UVRC.MCC
UVRC.MCC.MVinput		= UVRC.MV_in;
UVRC.MCC.derivFilt 		= eps;
UVRC.MCC.AWBTintterm 	= -UVRC.MCC.MVinput; %UVRC.MCC.AWBTintterm 	= eps;
%{
UVRC.x0 = [...
	UVRC.scrb(1).P; UVRC.scrb(1).T; ...
	UVRC.scrb(2).P; UVRC.scrb(2).T; ...
	UVRC.asc(1).MVinput; UVRC.asc(2).MVinput; ...
	UVRC.TIC(1).AWBTintterm; UVRC.TIC(1).derivFilt; UVRC.TIC(1).MVinput; ...
	UVRC.TIC(2).AWBTintterm; UVRC.TIC(2).derivFilt; UVRC.TIC(2).MVinput; ...
	UVRC.MCC.AWBTintterm; UVRC.MCC.derivFilt; UVRC.MCC.MVinput ...
	];
%}
UVRC.x0 = [...
	UVRC.scrb(1).P; UVRC.scrb(1).T; ...
	UVRC.scrb(2).P; UVRC.scrb(2).T; ...
	UVRC.TIC(1).AWBTintterm;  ...
	UVRC.TIC(2).AWBTintterm; ...
	UVRC.MCC.AWBTintterm ...
	];

%PI_pos LIC_001
LIC_001.MVinput			= LV_001.c;
LIC_001.derivFilt 		= eps;
LIC_001.AWBTintterm 	= -LIC_001.MVinput; %LIC_001.AWBTintterm 	= eps;
%LIC_001.x0 = [LIC_001.AWBTintterm; LIC_001.derivFilt; LIC_001.MVinput];
LIC_001.x0 = [LIC_001.AWBTintterm];

%PI_pos LIC_002
LIC_002.MVinput			= LV_002.c;
LIC_002.derivFilt 		= eps;
LIC_002.AWBTintterm 	= -LIC_002.MVinput; %LIC_002.AWBTintterm 	= eps;
%LIC_002.x0 = [LIC_002.AWBTintterm; LIC_002.derivFilt; LIC_002.MVinput];
LIC_002.x0 = [LIC_002.AWBTintterm];

%PI_pos LIC_003
LIC_003.MVinput			= LV_003.c;
LIC_003.derivFilt 		= eps;
LIC_003.AWBTintterm 	= -LIC_003.MVinput; %LIC_003.AWBTintterm 	= eps;
%LIC_003.x0 = [LIC_003.AWBTintterm; LIC_003.derivFilt; LIC_003.MVinput];
LIC_003.x0 = [LIC_003.AWBTintterm];

%PI_pos LIC_004
LIC_004.MVinput			= LV_004.c;
LIC_004.derivFilt 		= eps;
LIC_004.AWBTintterm 	= -LIC_004.MVinput; %LIC_004.AWBTintterm 	= eps;
%LIC_004.x0 = [LIC_004.AWBTintterm; LIC_004.derivFilt; LIC_004.MVinput];
LIC_004.x0 = [LIC_004.AWBTintterm];

%PI_pos LIC_005
LIC_005.MVinput			= LV_005.c;
LIC_005.derivFilt 		= eps;
LIC_005.AWBTintterm 	= -LIC_005.MVinput; %LIC_005.AWBTintterm 	= eps;
%LIC_005.x0 = [LIC_005.AWBTintterm; LIC_005.derivFilt; LIC_005.MVinput];
LIC_005.x0 = [LIC_005.AWBTintterm];

%PI_pos LIC_006
LIC_006.MVinput			= LV_006.c;
LIC_006.derivFilt 		= eps;
LIC_006.AWBTintterm 	= -LIC_006.MVinput; %LIC_006.AWBTintterm 	= eps;
%LIC_006.x0 = [LIC_006.AWBTintterm; LIC_006.derivFilt; LIC_006.MVinput];
LIC_006.x0 = [LIC_006.AWBTintterm];

%PI_pos PIC_001
PIC_001.MVinput			= PV_001.c;
PIC_001.derivFilt 		= eps;
PIC_001.AWBTintterm 	= -PIC_001.MVinput; %PIC_001.AWBTintterm 	= eps;
%PIC_001.x0 = [PIC_001.AWBTintterm; PIC_001.derivFilt; PIC_001.MVinput];
PIC_001.x0 = [PIC_001.AWBTintterm];

%PI_pos TIC_001
TIC_001.MVinput			= (P_002.Q_in-5e3)/(25e3-5e3);
TIC_001.derivFilt 		= eps;
TIC_001.AWBTintterm 	= TIC_001.MVinput; %TIC_001.AWBTintterm 	= eps;
%TIC_001.x0 = [TIC_001.AWBTintterm; TIC_001.derivFilt; TIC_001.MVinput];
TIC_001.x0 = [TIC_001.AWBTintterm];

%PI_pos TIC_002
TIC_002.MVinput			= P_003.Q_out/2.5e3;
TIC_002.derivFilt 		= eps;
TIC_002.AWBTintterm 	= -TIC_002.MVinput; %TIC_002.AWBTintterm 	= eps;
%TIC_002.x0 = [TIC_002.AWBTintterm; TIC_002.derivFilt; TIC_002.MVinput];
TIC_002.x0 = [TIC_002.AWBTintterm];

%----- Steady-state computation:
x0 = [...
	SG_001.x0; ...
	V_TO_001.x0; ...
	TO_001.x0; ...
	V_TO_002.x0; ...
	TO_002.x0; ...
	UVRC.x0; ...
	V_001.x0; ...
	LIC_001.x0; ...
	LIC_002.x0; ...
	LIC_003.x0; ...
	LIC_004.x0; ...
	LIC_005.x0; ...
	LIC_006.x0; ...
	PIC_001.x0; ...
	TIC_001.x0; ...
	TIC_002.x0  ...
	];

u0 = [...
	feed.Qf_in; feed.Rgo; feed.BSW_in; feed.T_in; ...
	LV_002.P_out;
	PV_001.P_out;
	stk.P; ...
	P_001.DP_C; P_001.DP_H; ...
	P_002.DP; ...
	dslt.Qw; dslt.T; ...
	UVRC.cool(1).DP; ...
	UVRC.cool(2).DP; ...
	P_003.DP ...
	];

yss=[];

[x0,u0]=uotr_ctrl_ss(x0,u0,yss,par);
x0=[x0;eps];	%Adding state Vg,flare as null in t0
%}
%----- State-space linear model:
%sys=sfunc2ss(@uotr_ctrl,x0,u0,par)


%----- Dynamic simulation:
load('uotr_ctrl_ss.mat');
%tstart=0; tend=24*3600; tstep=tend; tspan=tstart:tstep:tend;
tstart=0; tend=7200; tstep=1; tspan=tstart:tstep:tend;
%tstart=0; tend=900; tstep=1; tspan=tstart:tstep:tend;
u=u0;               %Successful
%u(1)=0.999*u0(1);  %Successful
%u(1)=0.985*u0(1); 	%Successful
%u(1)=0.95*u0(1);   %Successful
%u(1)=0.9*u0(1);    %Successful
%u(1)=0.85*u0(1);   %Successful
%u(1)=0.80*u0(1);   %Successful
%u(1)=0.70*u0(1);   %Successful
%u(1)=0.60*u0(1);   %Successful
u(1)=0.50*u0(1);   %Successful
%u(1)=0.40*u0(1);    %Fail
%u(1)=0.30*u0(1);
%u(1)=0.20*u0(1);
%u(1)=0.10*u0(1);   %Fail

%Jacobian matrix
Nx=size(x0,1);
fhand=@(t,x)uotr_ctrl(t,x,u0,1,[],par);
fval=fhand(0,x0);
[dFdx,~]=numjac(fhand,0,x0,fval,zeros(Nx,1));

Nu=size(u0,1);
fhand=@(t,u)uotr_ctrl(t,x0,u,1,[],par);
fval=fhand(0,u0);
[dFdu,~]=numjac(fhand,0,u0,fval,zeros(Nu,1));

%opt=odeset('RelTol',1e-5,'AbsTol',1e-8,'Stats','on','Jacobian',dFdx);
opt=odeset('RelTol',1e-5,'AbsTol',1e-8,'Stats','on');

tic;
[t,x]=ode15s(@(t,x)uotr_ctrl(t,x,u,1,[],par),tspan,x0,opt);
%[~,xlin]=ode15s(@(t,x)uotr_ctrl_lin(t,x,u,dFdx,dFdu,x0,u0),tspan,x0,opt);
xlin=modelo_eq(t,u(1)*1e-3,[t(1);t(end)],[]*1e-3,x0,u0,dFdx,dFdu);
timer=toc;
figure(1); plot(t,x(:,6),'b-',t,x(:,10),'r-',t,x(:,16),'g-',t,xlin(:,6),'b--',t,xlin(:,10),'r--',t,xlin(:,16),'g--'); %[SG_001.P; V_TO_001.P; V_TO_002.P]
figure(2); plot(t,x(:,end),'b-',t,xlin(:,end),'b--');
%figure(2); plot(t,x(:,[7;11;17;27]),'-',t,xlin(:,[7;11;17;27]),'--'); %[SG_001.T; V_TO_001.T; V_TO_002.T, V_002.T]
%figure(3); plot(t,x(:,[4;8;14]),'-',t,xlin(:,[4;8;14]),'--'); %[SG_001.hCL; V_TO_001.hL; V_TO_002.hL]
%figure(4); plot(t,x(:,[1;12;18]),'-',t,xlin(:,[1;12;18]),'--'); %[SG_001.hWCS; TO_001.hW; V_TO_002.hW]
%figure(5); plot(t,x(:,[24;25])); %[UVRC.asc(1).MVinput; UVRC.asc(2).MVinput]
%}