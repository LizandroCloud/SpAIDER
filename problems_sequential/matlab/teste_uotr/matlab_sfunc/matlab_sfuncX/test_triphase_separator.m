%Script to run triphase_separator Level-1 MATLAB S-Function
clear

%Physical properties
global PP
PP=blackoil;

%Setting model parameters
par=cell(8,1);
par{1}=141.5/(30+131.5);          %d60, Oil specific gravity []
par{2}=(1/1.8)*(11.8*par{1})^3;   %Tbn, Oil normal boiling point [K]
par{3}=255;     				  %Cws, Water salt concentration [g/L]
par{4}=0.7;                       %dg, Gas specific gravity []
par{5}=16.0;                      %LCS, Separation chamber (CS) length [m]
par{6}=4.00;                      %LCL, Lights chamber (CL) length [m]
par{7}=5.00;                      %D, Vessel diameter [m]
par{8}=2.50;                      %hweir, Weir height [m]

%Model states linearization point
hWCS   	=   0.5; 			%Water phase level in CS          	[m]
BSWLCS 	=   0.005;			%Basic sediments & water in CS 		[]
hCS    	=   1.001*par{8};	%Basic sediments & water in CS 		[]
hCL    	=   1.0;			%Liquid level in CL                 [m]
BSWLCL	=	0.005;			%Basic sediments & water in CL 		[]
P       =   2101.325;		%Pressure                           [kPa]
T       =   300;			%Temperature                        [K]

%Model inputs linearization point
Qf_in=24000;              	%Inlet liquid std. vol. flow rate  [Sm^3/d]
BSW_in=0.2;               	%Inlet basic sediments & water     []
Qg_in=200*Qf_in*(1-BSW_in); %Inlet gas std. vol. flow rate     [Sm^3/d]
T_in=300;                   %Inlet temperature                 [K]
effWL=0.4;
Qf=Qf_in*(1-BSW_in*effWL);  	%Liquid std. vol. outlet flow rate   [Sm^3/d]
Qw=Qf_in*BSW_in*effWL;        	%Water std. vol. outlet flow rate    [Sm^3/d]
Qg=Qg_in;                 		%Gas std. vol. outlet flow rate      [Sm^3/d]

%Getting state-space linear model
cinit=[hWCS; BSWLCS; hCS; hCL; BSWLCL; P; T];
u0=[Qf_in; Qg_in; BSW_in; T_in; Qf; Qw; Qg];
sys=sfunc2ss(@triphase_separator,cinit,u0,par)