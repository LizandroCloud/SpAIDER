%Script to run uvrc_ctrl Level-1 MATLAB S-Function
clear, close all

%Physical properties
global PP
PP=blackoil;
%{
%----- Setting model parameters:
%Global UVRC parameters
par.d60		= [141.5/(30+131.5), 141.5/(30+131.5)]; %Oil specific gravity []
par.dg 		= [1.05,1.5];	%Gas specific gravity 				[]
par.Nr_Max 	= 12000;		%Maximum shitf assemble rotation 	[rpm]
par.nS 		= 2;			%Number of compression stages 		[]
par.Qout_Max= [2e3, 2e3];	%Max heat removed from system 		[kW]

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
feed(1).Qg 		= 1.5e5;				%Pressure                   			[kPa]
feed(2).Qg 		= 1.5e5;				%Pressure                   			[kPa]
feed(1).T 		= 360;					%Temperature                   			[K]
feed(2).T 		= 360;					%Temperature                   			[K]
cool(1).Q_out 	= 0.4e3;				%Heat removed from system 				[kW]
cool(2).Q_out 	= 1.25e3;				%Heat removed from system 				[kW]
cool(1).DP 		= 30;					%Pressure drop 							[kPa]
cool(2).DP 		= 30;					%Pressure drop 							[kPa]
P_out 			= 2101.325;				%Outlet pressure						[kPa]

%----- States linearization point:
%scrubber_gas uvrc.scrb(k)
scrb(1).P       = 301.325 - cool(1).DP; %Pressure                   			[kPa]
scrb(1).T       = 320;                  %Temperature                   			[K]
scrb(2).P       = 801.325 - cool(2).DP; %Pressure                   			[kPa]
scrb(2).T       = 320;                  %Temperature                   			[K]

scrb(1).x0 = [scrb(1).P; scrb(1).T];
scrb(2).x0 = [scrb(2).P; scrb(2).T];

%P_pos uvrc.asc(k)
asc(1).MVinput 		= 0.2;
asc(2).MVinput 		= 0.2;

asc(1).x0 = asc(1).MVinput;
asc(2).x0 = asc(2).MVinput;

%PID_pos uvrc.TIC(k)
TIC(1).MVinput			= cool(1).Q_out/par.Qout_Max(1);
TIC(1).derivFilt 		= eps;
TIC(1).AWBTintterm 		= eps;%TIC(1).AWBTintterm 		= TIC(1).MVinput;

TIC(2).MVinput			= cool(2).Q_out/par.Qout_Max(2);
TIC(2).derivFilt 		= eps;
TIC(2).AWBTintterm 		= eps;%TIC(2).AWBTintterm 		= TIC(2).MVinput;

TIC(1).x0 = [TIC(1).AWBTintterm; TIC(1).derivFilt; TIC(1).MVinput];
TIC(2).x0 = [TIC(2).AWBTintterm; TIC(2).derivFilt; TIC(2).MVinput];

%PID_pos uvrc.MCC
MCC.MVinput				= 0.95;
MCC.derivFilt 			= eps;
MCC.AWBTintterm 		= eps;%MCC.AWBTintterm 		= MCC.MVinput;

MCC.x0 = [MCC.AWBTintterm; MCC.derivFilt; MCC.MVinput];

%----- Steady-state computation:
x0 = [...
	scrb(1).x0; scrb(2).x0; ...
	asc(1).x0; asc(2).x0; ...
	TIC(1).x0; TIC(2).x0; ...
	MCC.x0];

u0 = [...
	feed(1).Qg; feed(1).T; ...
	feed(2).Qg; feed(2).T; ...
	cool(1).DP; ...
	cool(2).DP; ....
	P_out];

yss = [];

[x0,u0]=uvrc_ctrl_ss(x0,u0,yss,par);
%}
%----- State-space linear model:
%{
for k=1:par.nS
u0(2*(k-1)+1)                   = feed(k).P;
u0(2*(k-1)+2)                   = feed(k).T;
u0(2*par.nS+(k-1)+1)            = cool(k).DP;
end
u0(2*par.nS+par.nS+par.nS+1)  	= P_out;
u0 = u0(:);
%}
%sys=sfunc2ss(@uvrc_ctrl,x0,u0,par)


%----- Dynamic simulation:
load('uvrc_ctrl_ss.mat');
tstart=0; tend=7200; tstep=1; tspan=tstart:tstep:tend;
u=u0;
u(end)=100+u0(end);

tic;
opt=odeset('RelTol',1e-5,'AbsTol',1e-8,'Stats','on');
[t,x]=ode15s(@(t,x)uvrc_ctrl(t,x,u,1,[],par),tspan,x0,opt);
timer=toc;
figure(1); plot(t,x(:,[1;3])); %[scrb(1).P; scrb(2).P]
figure(2); plot(t,x(:,[2;4])); %[scrb(1).T; scrb(2).T]
figure(3); plot(t,x(:,[5;6])); %[asc(1).MVinput; asc(2).MVinput]
figure(4); plot(t,x(:,[9;12])); %[TIC(1).MVinput; TIC(2).MVinput]
figure(5); plot(t,x(:,end)*par.Nr_Max); %[comp.Nr]
%}