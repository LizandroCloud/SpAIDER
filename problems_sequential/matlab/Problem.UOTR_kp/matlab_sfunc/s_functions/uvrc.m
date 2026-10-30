function [sys,x0,str,ts,simStateCompliance]=uvrc(t,x,u,flag,cinit,par)
%uvrc - Vapor recovery unit non-linear model S-Function
%
%   This model describes the dynamic behavior of an upstream topside
%	facility vapor recovery unit, using the black-oil approach.
%
%   [sys,x0,str,ts,simStateCompliance]=uvrc(t,x,u,flag,cinit,par)
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
%par.scrb(k): 			{C,D} 					k-th gas scrubber model parameters
%par.comp(k): 			{Fgm,effM,PrCC,EffSPOL} k-th compressor model parameters
%par.ascv(k): 			{Cv,Fgm,DPrT}    		k-th anti-surge control valve parameters
%--------------------------------------------------------------------
%----- States -------------------------------------------------------
%x(2*(k-1)+1):  		scrb(k).P     			Pressure       						[kPa]
%x(2*(k-1)+2):  		scrb(k).T     			Temperature       					[K]
%--------------------------------------------------------------------
%----- Inputs -------------------------------------------------------
%u(2*(k-1)+1):  		feed(k).Qg     			Feed gas std. vol. flow rate    	[Sm^3/d]
%u(2*(k-1)+2):  		feed(k).T     			Feed temperature       				[K]
%u(2*nS+2*(k-1)+1):		cool(k).Q_out   		Heat removed from the system        [kW]
%u(2*nS+2*(k-1)+2):		cool(k).DP   			Pressure drop        				[kPa]
%u(2*nS+2*nS+1):		P_out   				Outlet pressure        				[kPa]
%u(2*nS+2*nS+2):		MV_in     				Input signal from controller        []
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

%   2023. Version 0. Coded by Diego Trica (diegotrica@gmail.com)

%Blackoil object class property package. On the main script the
%user must declare PP as global and PP=blackoil.
global PP

%----- Model parameters:
d60         = par.d60;   	%Oil specific gravity []
dg          = par.dg;   	%Gas specific gravity []
Nr_Max      = par.Nr_Max;   %Maximum shift rotation [rpm]
nS          = par.nS;   	%Number of compression stages

for k=1:nS
mix(k).nG = 3;
mix(k).par	= {dg(k),mix(k).nG};
cool(k).par	= {dg(k)};
scrb(k).par	= {d60(k),dg(k),par.scrb(k).C,par.scrb(k).D};
comp(k).par = {dg(k),par.comp(k).Fgm,Nr_Max,par.comp(k).PrCC,par.comp(k).EffSPOL};
FIT(k).par	= {0,1e6};		%{RangeMin,RangeMax} Sm^3/d
asc(k).par 	= {'Reverse'};	%{action}
ascv(k).par = {dg(k),par.ascv(k).Cv,par.ascv(k).Fgm,par.ascv(k).DPrT};
end

%----- States:	
for k=1:nS
%scrubber_gas scrb(k)
scrb(k).x 		= x(2*(k-1)+1:2*(k-1)+2);
scrb(k).P     	= scrb(k).x(1);
scrb(k).T     	= scrb(k).x(2);
end

%----- Inputs:
for k=1:nS
feed(k).Qg     	= u(2*(k-1)+1);
feed(k).T     	= u(2*(k-1)+2);
cool(k).Q_out	= u(2*nS+2*(k-1)+1);
cool(k).DP		= u(2*nS+2*(k-1)+2);
asc(k).mode		= 0;	%auto
asc(k).KC		= 2.0; %1.0;
asc(k).TI		= 30.0;%90.0;	%[s]
asc(k).TD		= 0.0;	%[s] (PI controller)
end
P_out 			= u(2*nS+2*nS+1);
MV_in 			= u(2*nS+2*nS+2);

%----- Internal calculations:
%Stage discharge pressure [kPa]
stg(1:nS-1).Pd = scrb(2:nS).P + cool(2:nS).DP;
stg(nS).Pd = P_out;

for k=1:nS
%compressor comp(k)
comp(k).P_in 	= scrb(k).P;
comp(k).T_in 	= scrb(k).T;
comp(k).P_out 	= stg(k).Pd;
comp(k).Nr		= MV_in*Nr_Max;
comp(k).u = [...
	comp(k).P_in; ...
	comp(k).T_in; ...
	comp(k).P_out; ...
	MV_in];
comp(k).y = compressor(t,[],comp(k).u,3,[],comp(k).par);
comp(k).Qg 		= comp(k).y(1);
comp(k).T 		= comp(k).y(2);
comp(k).W 		= comp(k).y(3);

%Anti-surge controller: Flow transmitter with dynamic set-point
%Flow transmitter
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
asc(k).y = P_pos(t,[],asc(k).u,3,[],asc(k).par);
asc(k).MVoutput = asc(k).y(1);

%valv_gas ascv(k)
ascv(k).P_in 	= stg(k).Pd;
ascv(k).T_in	= comp(k).T;
ascv(k).P_out	= scrb(k).P + cool(k).DP;
ascv(k).c		= asc(k).MVoutput;
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
	ascv(k).Qg = asc(k).MVoutput*comp(k).Qg;
	if k==1
		stg(k).Qg_th = 0;
		stg(k).T_th  = 300;
	else
		stg(k).Qg_th = (1-asc(k-1).MVoutput)*comp(k-1).Qg;
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

%cooler_gas cool(k)
cool(k).Qg_in 	= mix(k).Qg;
cool(k).P_in 	= scrb(k).P + cool(k).DP;
cool(k).T_in 	= mix(k).T;
cool(k).P_out 	= scrb(k).P;
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
sys(1)=2*nS;%Number of continuous states
sys(2)=0;%Number of discrete states
sys(3)=nS+3*nS+nS+2*nS;%Number of outputs
sys(4)=2*nS+2*nS+2;%Number of inputs
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
end

%dxdt allocation:
dxdt=[];
for k=1:nS
dxdt=[dxdt; ...
    scrb(k).dxdt];
end

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