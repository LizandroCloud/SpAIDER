%Script to run uotr Level-1 MATLAB S-Function
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
UVRC.Nr_Max 		= 12000;	%Maximum shitf assemble rotation 	[rpm]
UVRC.nS 			= 2;		%Number of compression stages 		[]

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
par.UVRC.nS		= UVRC.nS;
par.UVRC.scrb   = UVRC.scrb;
par.UVRC.comp   = UVRC.comp;
par.UVRC.ascv   = UVRC.ascv;
par.V_001 		= {V_001.C, V_001.D};

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
	
UVRC.asc(1).MVinput 	= 0.01;
UVRC.asc(2).MVinput 	= 0.1;

UVRC.x0 = [...
	UVRC.scrb(1).P; UVRC.scrb(1).T; ...
	UVRC.scrb(2).P; UVRC.scrb(2).T; ...
	UVRC.asc(1).MVinput; ...
	UVRC.asc(2).MVinput];
	
%triphase_separator SG_001 feed stream
feed.Qf_in		= 24000;    %Inlet liquid std. vol. flow rate  	[Sm^3/d]
feed.BSW_in		= 0.1;      %Inlet basic sediments & water     	[]
feed.Rgo 		= 250;		%Gas-oil ratio 						[Sm^3/Sm^3]
							%Inlet gas std. vol. flow rate 		[Sm^3/d]
feed.T_in		= 290;		%Inlet temperature                 	[K]

%valv_liq LV_001
LV_001.c		= 0.5;		%Input signal from controller      	[]	

%valv_water LV_002
LV_002.P_out	= 801.325;  %Outlet pressure                   	[kPa]
LV_002.c		= 0.5;		%Input signal from controller      	[]

%valv_gas PV_001
PV_001.P_out	= 501.325;	%Outlet pressure                   	[kPa]
PV_001.c		= 0.5;		%Input signal from controller      	[]

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
%{
%valv_gas PV_002
PV_002.P_out 	= 101.325;	%Outlet pressure                   	[kPa]
PV_002.c     	= 1e-5;		%Input signal from controller      	[]
%}
%mixer_W3ph D_001 desalt water stream
dslt.Qw 		= 1200;		%Std. vol. flow rate 				[Sm^3/d]
dslt.T 			= 360;		%Temperature 						[K]

%valv_liq LV_005
LV_005.c     	= 0.5;		%Input signal from controller      	[]

%valv_water LV_006
LV_006.c     	= 0.5;		%Input signal from controller      	[]

%Stock tank stream
stk.P = 101.325;			%Stock tank pressure 				[kPa]
%{
%valv_gas PV_003
PV_003.P_out 	= 101.325;	%Outlet pressure                   	[kPa]
PV_003.c     	= 1e-5;		%Input signal from controller      	[]
%}
%vapor recovery unit UVRC
UVRC.MV_in 		= 0.85;		%Input signal from capacity control	[]

%----- Steady-state computation:
x0 = [...
	SG_001.x0; ...
	V_TO_001.x0; ...
	TO_001.x0; ...
	V_TO_002.x0; ...
	TO_002.x0; ...
	UVRC.x0; ...
	V_001.x0];

u0 = [...
	feed.Qf_in; feed.Rgo; feed.BSW_in; feed.T_in; ...
	LV_001.c; ...
	LV_002.P_out; LV_002.c; ...
	PV_001.P_out; PV_001.c; ...
	stk.P; ...
	P_001.DP_C; P_001.DP_H; ...
	P_002.DP; P_002.Q_in; ...
	LV_003.c; ...
	LV_004.c; ...
	dslt.Qw; dslt.T; ...
	LV_005.c; ...
	LV_006.c; ...
	UVRC.cool(1).Q_out; UVRC.cool(1).DP; ...
	UVRC.cool(2).Q_out; UVRC.cool(2).DP; ...
	UVRC.MV_in; ...
	P_003.Q_out; P_003.DP];

yss=[];

[x0,u0]=uotr_ss(x0,u0,yss,par);
x0=[x0;eps];	%Adding state PV_001.Vg as null in t0
%}

%----- State-space linear model:
%sys=sfunc2ss(@uotr,x0,u0,par)


%----- Dynamic simulation:
load('uotr_ss.mat');
%tstart=0; tend=24*3600; tstep=tend; tspan=tstart:tstep:tend;
tstart=0; tend=7200; tstep=1; tspan=tstart:tstep:tend;
%tstart=0; tend=900; tstep=1; tspan=tstart:tstep:tend;
u=u0;
%u(1)=0.999*u0(1);                                  %Successful
%u(1)=0.985*u0(1); 									%Successful
%u(1)=0.95*u0(1); u(9)=0.9*u0(9); 					%Successful
%u(1)=0.9*u0(1);  u(7)=0.93*u0(7); u(9)=0.88*u0(9); %Successful
%u(1)=0.85*u0(1); u(7)=0.90*u0(7); u(9)=0.90*u0(9); %Fail
%u([1;5;7;15;19;9])=[0.88;0.97;0.95;0.98;0.95;0.88].*u0([1;5;7;15;19;9]); %Successful
%u([1;5;7;15;19;9])=[0.86;0.97;0.95;0.98;0.93;0.86].*u0([1;5;7;15;19;9]); %Successful
%u([1;5;7;15;19;9])=[0.82;0.96;0.95;0.98;0.95;0.80].*u0([1;5;7;15;19;9]); %Successful
%u([1;5;7;15;19;9])=[0.78;0.96;0.95;0.98;0.92;0.78].*u0([1;5;7;15;19;9]); %Successful	
%u([1;5;7;15;19;9])=[0.72;0.96;0.95;0.98;0.92;0.72].*u0([1;5;7;15;19;9]); %Successful	
%u([1;5;7;15;19;9])=[0.70;0.95;0.95;0.96;0.90;0.72].*u0([1;5;7;15;19;9]);  %Fail
%u([1;5;7;15;19;9])=[0.60;0.93;0.95;0.96;0.87;0.58].*u0([1;5;7;15;19;9]); %Successful	
%u([1;5;7;15;19;9;14;26])=[0.50;0.92;0.95;0.96;0.87;0.50;0.9;0.9].*u0([1;5;7;15;19;9;14;26]); %Fail

tic;
opt=odeset('RelTol',1e-5,'AbsTol',1e-8,'Stats','on');
[t,x]=ode15s(@(t,x)uotr(t,x,u,1,[],par),tspan,x0,opt);
timer=toc;
figure(1); plot(t,x(:,[6;10;16])); %[SG_001.P; V_TO_001.P; V_TO_002.P]
figure(2); plot(t,x(:,[7;11;17;27])); %[SG_001.T; V_TO_001.T; V_TO_002.T, V_002.T]
figure(3); plot(t,x(:,[4;8;14])); %[SG_001.hCL; V_TO_001.hL; V_TO_002.hL]
figure(4); plot(t,x(:,[1;12;18])); %[SG_001.hWCS; TO_001.hW; V_TO_002.hW]
figure(5); plot(t,x(:,[24;25])); %[UVRC.asc(1).MVinput; UVRC.asc(2).MVinput]
%}