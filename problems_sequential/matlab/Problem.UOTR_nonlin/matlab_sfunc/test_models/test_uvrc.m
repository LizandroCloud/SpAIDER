%Script to run uvrc Level-1 MATLAB S-Function
clear

%Physical properties
global PP
PP=blackoil;

%----- Setting model parameters:
%Global UVRC parameters
par.d60		= [141.5/(30+131.5), 141.5/(30+131.5)]; %Oil specific gravity []
par.dg 		= [1.05,1.5];	%Gas specific gravity 				[]
par.Nr_Max 	= 12000;		%Maximum shitf assemble rotation 	[rpm]
par.nS 		= 2;			%Number of compression stages 		[]

%compressor uvrc.scrb(k)
par.scrb(1).C 		= 3.0;				%Scrubber length 						[m]
par.scrb(1).D 		= 1.25;				%Scrubber diameter 						[m]

par.scrb(2).C 		= 3.0;				%Scrubber length 						[m]
par.scrb(2).D 		= 1.25;				%Scrubber diameter 						[m]

%compressor uvrc.comp(k)
par.comp(1).Fgm     = 1.2;     			%Specific heat capacity ratio factor 	[]
par.comp(1).PrCC	= [0.35,4.5,1.5];  	%Pres. ratio charac. curve par.   		[]
par.comp(1).EffSPOL	= [0.35,1.1,7.0];   %Isentropic eff. polynomial par.   		[]

par.comp(2).Fgm     = 1.2;     			%Specific heat capacity ratio factor 	[]
par.comp(2).PrCC	= [0.35,4.5,1.5];  	%Pres. ratio charac. curve par.   		[]
par.comp(2).EffSPOL	= [0.35,1.1,7.0];   %Isentropic eff. polynomial par.   		[]

%valv_gas uvrc.ascv(k)
par.ascv(1).Cv		= 300;              %Valve constant US units 				[gpm]
par.ascv(1).Fgm		= 1.2;              %Spec. heat cap. ratio factor 			[]
par.ascv(1).DPrT	= 0.70;     		%Pressure drop ratio factor 			[]
	
par.ascv(2).Cv		= 300;              %Valve constant US units 				[gpm]
par.ascv(2).Fgm		= 1.2;              %Spec. heat cap. ratio factor 			[]
par.ascv(2).DPrT	= 0.70;     		%Pressure drop ratio factor 			[]

%----- Inputs linearization point:
feed(1).Qg 		= 1.5e5;						%Pressure                   			[kPa]
feed(2).Qg 		= 1.5e5;						%Pressure                   			[kPa]
feed(1).T 		= 360;							%Temperature                   			[K]
feed(2).T 		= 360;							%Temperature                   			[K]
cool(1).Q_out 	= 0.4e3;						%Heat removed from system 				[kW]
cool(2).Q_out 	= 1.25e3;						%Heat removed from system 				[kW]
cool(1).DP 		= 30;							%Pressure drop 							[kPa]
cool(2).DP 		= 30;							%Pressure drop 							[kPa]
%ascv(1).c 		= 0.1;							%Valve stem position 					[]
%ascv(2).c 		= 0.1;							%Valve stem position 					[]
P_out 			= 2101.325;						%Outlet pressure						[kPa]
MV_in 			= 0.95;							%Input signal from capacity control		[]

%----- States linearization point:
scrb(1).P       = 301.325 - cool(1).DP;			%Pressure                   			[kPa]
scrb(1).T       = 320;                          %Temperature                   			[K]
scrb(2).P       = 801.325 - cool(2).DP;       	%Pressure                   			[kPa]
scrb(2).T       = 320;                          %Temperature                   			[K]

scrb(1).x0 = [scrb(1).P; scrb(1).T];
scrb(2).x0 = [scrb(2).P; scrb(2).T];

asc(1).AWBTintterm 	= 0;
asc(1).derivFilt 	= 0;
asc(1).MVinput 		= 0.1;
asc(2).AWBTintterm 	= 0;
asc(2).derivFilt 	= 0;
asc(2).MVinput 		= 0.01;

asc(1).x0 = [asc(1).AWBTintterm; asc(1).derivFilt; asc(1).MVinput];
asc(2).x0 = [asc(2).AWBTintterm; asc(2).derivFilt; asc(2).MVinput];

%----- Steady-state computation:
x0 = [...
	scrb(1).x0; scrb(2).x0; ...
	asc(1).x0; asc(2).x0];

u0 = [...
	feed(1).Qg; feed(1).T; ...
	feed(2).Qg; feed(2).T; ...
	cool(1).Q_out; cool(1).DP; ...
	cool(2).Q_out; cool(2).DP; ....
	%ascv(1).c; ascv(2).c; ...
	P_out; MV_in];

yss = [];

[x0,u0]=uvrc_ss(x0,u0,yss,par);

%----- State-space linear model:,
%{
for k=1:par.nS
u0(2*(k-1)+1)                   = feed(k).P;
u0(2*(k-1)+2)                   = feed(k).T;
u0(2*par.nS+2*(k-1)+1)          = cool(k).Q_out;
u0(2*par.nS+2*(k-1)+2)          = cool(k).DP;
u0(2*par.nS+2*par.nS+(k-1)+1)   = ascv(k).c;
end
u0(2*par.nS+2*par.nS+par.nS+1)  = P_out;
u0(2*par.nS+2*par.nS+par.nS+2)  = MV_in;
u0 = u0(:);
%}
%sys=sfunc2ss(@uvrc,x0,u0,par)