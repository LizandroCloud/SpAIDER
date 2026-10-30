function [sys,x0,str,ts,simStateCompliance]=uvrc_ctrl(t,x,u,flag,cinit,par)
%uvrc - Vapor recovery unit non-linear model S-Function with PID 
%controllers
%
%   This model describes the dynamic behavior of an upstream topside
%	facility vapor recovery unit, using the black-oil approach, with
%	PID controllers for process variables such as levels, pressures
%	and temperatures.
%
%   [sys,x0,str,ts,simStateCompliance]=uvrc_ctrl(t,x,u,flag,cinit,par)
%
%   sys is the S-function output, which depends on which flag is 
%   being inputted by simulink. x0 are the states initial condi-
%   tions. See the S-function Level-1 MATLAB reference page for more
%   details about str, ts and simStateCompliance S-Function ouputs.
%
%   t is the independent variable time. x are the state variables.
%   u are the input variables, combined into one vector of manipula-
%   ted and disturbs. flag is the Simulink variable inputted to get a
%   specific output from the S-function. cinit are variable set for 
%	state initial condition computation. par are the model parameters.
%
%----- Model parameters ---------------------------------------------
%par.d60: 				d60      				Oil specific gravity                []
%par.dg: 				dg      				Gas specific gravity                []
%par.Nr_Max: 			Nr_Max  				Maximum shitf assemble rotation     [rpm]
%par.nS: 				nS  					Number of compression stages 		[]
%par.Qout_Max(k): 		Qout_Max(k)  			Max heat removed from system 		[kW]
%par.scrb(k): 			{C,D} 					k-th gas scrubber model parameters
%par.comp(k): 			{Fgm,effM,PrCC,EffSPOL} k-th compressor model parameters
%par.ascv(k): 			{Cv,Fgm,DPrT}    		k-th anti-surge control valve parameters
%--------------------------------------------------------------------
%----- States -------------------------------------------------------
%x(2*(k-1)+1):  		scrb(k).P     			Pressure       						[kPa]
%x(2*(k-1)+2):  		scrb(k).T     			Temperature       					[K]
%x(2*nS+(k-1)+1):  		asc(k).MVinput     		MV tracking       					[]
%x(2*nS+nS+3*(k-1)+1):  TIC(k).AWBTintterm 		Integral action with AWBT			[]
%x(2*nS+nS+3*(k-1)+2): 	TIC(k).derivFilt   		Derivative action with filter      	[]
%x(2*nS+nS+3*(k-1)+3): 	TIC(k).MVinput     		MV tracking       					[]
%x(2*nS+nS+3*nS+1):  	MCC.AWBTintterm 		Integral action with AWBT			[]
%x(2*nS+nS+3*nS+2):  	MCC.derivFilt   		Derivative action with filter      	[]
%x(2*nS+nS+3*nS+3):  	MCC.MVinput     		MV tracking       					[]
%--------------------------------------------------------------------
%----- Inputs -------------------------------------------------------
%u(2*(k-1)+1):  		feed(k).Qg     			Feed gas std. vol. flow rate    	[Sm^3/d]
%u(2*(k-1)+2):  		feed(k).T     			Feed temperature       				[K]
%u(2*nS+(k-1)+1):		cool(k).DP   			Pressure drop        				[kPa]
%u(2*nS+nS+1):			P_out   				Outlet pressure        				[kPa]
%--------------------------------------------------------------------
%------ Outputs -----------------------------------------------------
%y((k-1)+1): 			feed(k).P      			Feed pressure    					[kPa]
%y(nS+3*(k-1)+1): 		comp(k).Qg      		Gas std. vol. flow rate             [Sm^3/d]
%y(ns+3*(k-1)+2): 		comp(k).T       		Gas std. vol. flow rate             [Sm^3/d]
%y(ns+3*(k-1)+3): 		comp(k).W				Compressor power 				 	[kW]
%y(ns+3*nS+(k-1)): 		stg(k).Qg_th			Gas std. vol. flow rate             [Sm^3/d]
%y(nS+3*nS+nS+2*(k-1)+1): scrb(k).P 			Pressure                         	[kPa]
%y(nS+3*nS+nS+2*(k-1)+2): scrb(k).T 			Temperature                         [K]
%--------------------------------------------------------------------
%------ Initial condition set ---------------------------------------
%cinit = x
%--------------------------------------------------------------------

%   2024-06-09. Version 0. Coded by Diego Trica (diegotrica@gmail.com)

%Blackoil object class property package. On the main script the
%user must declare PP as global and PP=blackoil.
global PP

%----- Model parameters:
d60         	= par.d60;   			%Oil specific gravity 			[]
dg          	= par.dg;   			%Gas specific gravity 			[]
Nr_Max      	= par.Nr_Max;   		%Maximum shift rotation 		[rpm]
nS          	= par.nS;   			%Number of compression stages
Qout_Max(1:nS) 	= par.Qout_Max(1:nS);	%Max heat removed from system 	[kW]

for k=1:nS
mix(k).nG = 3;
mix(k).par	= {dg(k),mix(k).nG};
cool(k).par	= {dg(k)};
scrb(k).par	= {d60(k),dg(k),par.scrb(k).C,par.scrb(k).D};
comp(k).par = {dg(k),par.comp(k).Fgm,Nr_Max,par.comp(k).PrCC,par.comp(k).EffSPOL};
FIT(k).par	= {0,1e6};					%{RangeMin,RangeMax} [Sm^3/d]
asc(k).par 	= {'Reverse'};				%{action}
ascv(k).par = {dg(k),par.ascv(k).Cv,par.ascv(k).Fgm,par.ascv(k).DPrT};
TIT(k).par	= {30+273.15,300+273.15}; 	%{RangeMin,RangeMax} [K]
TIC(k).par 	= {'Direct'};				%{action}
end
PIT.par		= {0,500};					%{RangeMin,RangeMax} [kPa]
MCC.par 	= {'Direct'};				%{action}

%----- States:	
for k=1:nS
%scrubber_gas scrb(k)
scrb(k).x 		= x(2*(k-1)+1:2*(k-1)+2);
scrb(k).P     	= scrb(k).x(1);
scrb(k).T     	= scrb(k).x(2);

%P_pos asc(k)
asc(k).x 		= x(2*nS+(k-1)+1);
asc(k).MVinput  = asc(k).x(1);

%PID_pos TIC(k)
TIC(k).x 			= x(2*nS+nS+3*(k-1)+1:2*nS+nS+3*(k-1)+3);
TIC(k).AWBTintterm 	= TIC(k).x(1);
TIC(k).derivFilt   	= TIC(k).x(2);
TIC(k).MVinput     	= TIC(k).x(3);
end

%PID_pos MCC
MCC.x 				= x(2*nS+nS+3*nS+1:2*nS+nS+3*nS+3);
MCC.AWBTintterm 	= MCC.x(1);
MCC.derivFilt   	= MCC.x(2);
MCC.MVinput     	= MCC.x(3);

%----- Inputs:
for k=1:nS
feed(k).Qg     	= u(2*(k-1)+1);
feed(k).T     	= u(2*(k-1)+2);
cool(k).DP		= u(2*nS+(k-1)+1);
end
P_out 			= u(2*nS+nS+1);

%----- Control parameters:
for k=1:nS
%asc(k) anti-surge control
asc(k).mode		= 0;	%auto
asc(k).KC		= 2.0;
asc(k).TI		= 30.0;	%[s]
asc(k).TD		= 0.0;	%[s] (PI controller)

%cool(k) temperature control
TIC(k).mode		= 0;	%auto
TIT(k).SPinput 	= 320;	%[K]
TIC(k).KC		= 1.0;
TIC(k).TI		= 60.0;	%[s]
TIC(k).TD		= 0.0;	%[s] (PI controller)
end

%MCC: Master compressor capacity control
MCC.mode		= 0;	%auto
PIT.SPinput 	= 301.325-cool(1).DP;	%[kPa]
MCC.KC			= 0.5;
MCC.TI			= 600.0;%[s]
MCC.TD			= 0.0;	%[s] (PI controller)

%----- Internal calculations:
%Stage discharge pressure [kPa]
stg(1:nS-1).Pd = scrb(2:nS).P + cool(2:nS).DP;
stg(nS).Pd = P_out;

%MCC: Master compressor capacity control transmitter
PIT.PVinput = scrb(1).P;
PIT.u = [...
	PIT.PVinput; ...
	PIT.SPinput];
PIT.y = transmitter(t,[],PIT.u,3,[],PIT.par);
PIT.PVoutput = PIT.y(1);
PIT.SPoutput = PIT.y(2);

%PID_pos MCC
MCC.PV 			= PIT.PVoutput;
MCC.SP			= PIT.SPoutput;
MCC.u = [...
	MCC.mode; ...
	MCC.PV; ...
	MCC.SP; ...
	MCC.KC; MCC.TI; MCC.TD];
MCC.y = PID_pos(t,MCC.x,MCC.u,3,[],MCC.par);
MCC.MVoutput = MCC.y(1);

for k=1:nS
%compressor comp(k)
comp(k).P_in 	= scrb(k).P;
comp(k).T_in 	= scrb(k).T;
comp(k).P_out 	= stg(k).Pd;
comp(k).Nr		= MCC.MVoutput*Nr_Max; %comp(k).Nr		= MCC.MVinput*Nr_Max;
comp(k).u = [...
	comp(k).P_in; ...
	comp(k).T_in; ...
	comp(k).P_out; ...
	MCC.MVoutput];%MCC.MVinput];
comp(k).y = compressor(t,[],comp(k).u,3,[],comp(k).par);
comp(k).Qg 		= comp(k).y(1);
comp(k).T 		= comp(k).y(2);
comp(k).W 		= comp(k).y(3);

%Anti-surge controller: Flow transmitter with dynamic set-point
FIT(k).PVinput = comp(k).Qg;
if k==1
	FIT(k).SPinput = (dg(k)*PP.rho_AIRstd)\(3600*24)*... 
(2 + 2*0.5*(par.comp(k).PrCC(2)*1e-6)/(par.comp(k).PrCC(3)*1e-1)*comp(k).Nr*(comp(k).P_in/(100)));
end
if k==2
	FIT(k).SPinput = (dg(k)*PP.rho_AIRstd)\(3600*24)*... 
(3 + 3*0.5*(par.comp(k).PrCC(2)*1e-6)/(par.comp(k).PrCC(3)*1e-1)*comp(k).Nr*(comp(k).P_in/(100)));
end
FIT(k).u = [...
	FIT(k).PVinput; ...
	FIT(k).SPinput];
FIT(k).y = transmitter(t,[],FIT(k).u,3,[],FIT(k).par);
FIT(k).PVoutput = FIT(k).y(1);
FIT(k).SPoutput = FIT(k).y(2);

%P_pos asc(k)
asc(k).PV 		= FIT(k).PVoutput;
asc(k).SP		= FIT(k).SPoutput;
asc(k).u = [...
	asc(k).mode; ...
	asc(k).PV; ...
	asc(k).SP; ...
	asc(k).KC];	
asc(k).y = P_pos(t,asc(k).x,asc(k).u,3,[],asc(k).par);
asc(k).MVoutput = asc(k).y(1);

%valv_gas ascv(k)
ascv(k).P_in 	= stg(k).Pd;
ascv(k).T_in	= comp(k).T;
ascv(k).P_out	= scrb(k).P + cool(k).DP;
ascv(k).c		= asc(k).MVoutput;%ascv(k).c		= asc(k).MVinput;
ascv(k).u = [...
	ascv(k).P_in; ...
	ascv(k).T_in; ...
	ascv(k).P_out; ...
	1.00]; %ascv(k).c];
ascv(k).y = valv_gas(t,[],ascv(k).u,3,[],ascv(k).par);
%ascv(k).Qg 		= ascv(k).y(1);
ascv(k).T 		= ascv(k).y(2);
end

for k=1:nS
	ascv(k).Qg = asc(k).MVoutput*comp(k).Qg;%ascv(k).Qg = asc(k).MVinput*comp(k).Qg;
	if k==1
		stg(k).Qg_th = 0;
		stg(k).T_th  = 300;
	else
		stg(k).Qg_th = (1-asc(k-1).MVoutput)*comp(k-1).Qg;%stg(k).Qg_th = (1-asc(k-1).MVinput)*comp(k-1).Qg;
		stg(k).T_th  = comp(k-1).T;
	end
end

for k=1:nS
%mixer_G mix(k)
mix(k).G(1).Qg = feed(k).Qg;
mix(k).G(2).Qg = ascv(k).Qg;
mix(k).G(3).Qg = stg(k).Qg_th;

mix(k).G(1).T = feed(k).T;
mix(k).G(2).T = ascv(k).T;
mix(k).G(3).T = stg(k).T_th;

mix(k).u = [...
	mix(k).G(1).Qg; ...
	mix(k).G(1).T; ...
	mix(k).G(2).Qg; ...
	mix(k).G(2).T; ...
	mix(k).G(3).Qg; ...
	mix(k).G(3).T; ...
	(scrb(k).P + cool(k).DP)];
mix(k).y = mixer_G(t,[],mix(k).u,3,[],mix(k).par);
mix(k).Qg 		= mix(k).y(1);
mix(k).T 		= mix(k).y(2);

%cool(k) temperature transmitter
TIT(k).PVinput = scrb(k).T;
TIT(k).u = [...
	TIT(k).PVinput; ...
	TIT(k).SPinput];
TIT(k).y = transmitter(t,[],TIT(k).u,3,[],TIT(k).par);
TIT(k).PVoutput = TIT(k).y(1);
TIT(k).SPoutput = TIT(k).y(2);

%PID_pos TIC(k)
TIC(k).PV 		= TIT(k).PVoutput;
TIC(k).SP		= TIT(k).SPoutput;
TIC(k).u = [...
	TIC(k).mode; ...
	TIC(k).PV; ...
	TIC(k).SP; ...
	TIC(k).KC; TIC(k).TI; TIC(k).TD]; 
TIC(k).y = PID_pos(t,TIC(k).x,TIC(k).u,3,[],TIC(k).par);
TIC(k).MVoutput = TIC(k).y(1);

%cooler_gas cool(k)
cool(k).Qg_in 	= mix(k).Qg;
cool(k).P_in 	= scrb(k).P + cool(k).DP;
cool(k).T_in 	= mix(k).T;
cool(k).P_out 	= scrb(k).P;
cool(k).Q_out	= TIC(k).MVoutput*Qout_Max(k); %cool(k).Q_out	= TIC(k).MVinput*Qout_Max(k);
cool(k).u = [...
	cool(k).Qg_in; ...
	cool(k).P_in; ...
	cool(k).T_in; ...
	cool(k).P_out; ...
	cool(k).Q_out];
cool(k).y = cooler_gas(t,[],cool(k).u,3,[],cool(k).par);
cool(k).Qg 		= cool(k).y(1);
cool(k).T 		= cool(k).y(2);
end

switch flag,
%----- S-FUNCTION BUILDING PARAMETERS AND INITIAL CONDITIONS --------
    case 0,        
sys(1)=2*nS+nS+3*nS+3;
sys(2)=0;%Number of discrete states
sys(3)=nS+3*nS+nS+2*nS;%Number of outputs
sys(4)=2*nS+nS+1;%Number of inputs
sys(5)=0;%Reserved for root finding. Must be zero.
sys(6)=1;%Direct Feedthrough S-function parameter
sys(7)=1;%Number of sample times (at least one sample time is needed)

%----- States initial condition:
x0=cinit;
x0=x0(:);

%----- Other building parameters:
str=[];     %str is always an empty matrix
ts=[0 0];   %Sample time array [0 0] for continuous system
simStateCompliance='UnknownSimState';
%--------------------------------------------------------------------
%----- STATE DERIVATIVES --------------------------------------------
    case 1,
for k=1:nS
%scrubber_gas scrb(k)
scrb(k).Qg_in	= cool(k).Qg;
scrb(k).T_in	= cool(k).T;
scrb(k).Qg		= comp(k).Qg;
scrb(k).u = [...
    scrb(k).Qg_in; ...
    scrb(k).T_in; ...
    scrb(k).Qg];
scrb(k).dxdt = scrubber_gas(t,scrb(k).x,scrb(k).u,1,[],scrb(k).par);

%P_pos asc(k)
%{
asc(k).PV 		= FIT(k).PVoutput;
asc(k).SP		= FIT(k).SPoutput;
asc(k).u = [...
	asc(k).mode; ...
	asc(k).PV; ...
	asc(k).SP; ...
	asc(k).KC];	
%}	
asc(k).dxdt = P_pos(t,asc(k).x,asc(k).u,1,[],asc(k).par);

%PID_pos TIC(k)
%{
TIC(k).PV 		= TIT(k).PVoutput;
TIC(k).SP		= TIT(k).SPoutput;
TIC(k).u = [...
	TIC(k).mode; ...
	TIC(k).PV; ...
	TIC(k).SP; ...
	TIC(k).KC; TIC(k).TI; TIC(k).TD]; 
%}	
TIC(k).dxdt = PID_pos(t,TIC(k).x,TIC(k).u,1,[],TIC(k).par);
end

%PID_pos MCC
%{
MCC.PV 			= PIT.PVoutput;
MCC.SP			= PIT.SPoutput;
MCC.u = [...
	MCC.mode; ...
	MCC.PV; ...
	MCC.SP; ...
	MCC.KC; MCC.TI; MCC.TD]; 
%}
MCC.dxdt = PID_pos(t,MCC.x,MCC.u,1,[],MCC.par);

%dxdt allocation:
dxdt=[];
for k=1:nS
dxdt=[dxdt; ...
    scrb(k).dxdt];
end
for k=1:nS
dxdt=[dxdt; ...
	asc(k).dxdt];
end
for k=1:nS
dxdt=[dxdt; ...
	TIC(k).dxdt];
end
dxdt=[dxdt; MCC.dxdt];

%----- State derivatives for sys output flag=1
sys=dxdt;
sys=sys(:);
%--------------------------------------------------------------------
    case 2,
sys=[];
%----- OUPUTS -------------------------------------------------------
    case 3,
%y allocation
y=zeros(nS+3*nS+2*nS,1);
for k=1:nS
y((k-1)+1)					= cool(k).P_in;
y(nS+3*(k-1)+1)				= comp(k).Qg;
y(nS+3*(k-1)+2)				= comp(k).T;
y(nS+3*(k-1)+3)				= comp(k).W;
y(nS+3*nS+(k-1)+1)			= stg(k).Qg_th;
y(nS+3*nS+2*nS+2*(k-1)+1)	= scrb(k).P;
y(nS+3*nS+2*nS+2*(k-1)+2)	= scrb(k).T;
end

%----- Outputs for sys output flag=3:
sys=y;
sys=sys(:);
%--------------------------------------------------------------------
%----- S-FUNCTION OTHER FLAGS ---------------------------------------
    case 4,
sys=[];
    case 9,
sys=[];
    otherwise   % Unexpected flags
DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
end
%--------------------------------------------------------------------    
end