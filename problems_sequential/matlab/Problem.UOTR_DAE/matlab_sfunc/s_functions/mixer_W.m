function [sys,x0,str,ts,simStateCompliance]=mixer_W(t,x,u,flag,cinit,par)
%mixer_water - Water streams mixer non-linear model S-Function
%
%   This model describes the dynamic behavior of a mixer with water
%   brine only inlet streams for oil & gas upstream process, using 
%   the black-oil approach.
%
%   [sys,x0,str,ts,simStateCompliance]=mixer_W(t,x,u,flag,cinit,par)
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
%   Model parameters:
%   par{1}: Cws     Water salt concentration                	[g/L]
%   par{2}: nW      Number of W inlet streams                	[]
%
%   States:
%   There are no states for this model.
%
%   Inputs:
%	u(2*(k-1)+1): W(k).Qw  	Inlet W liquid std. vol. flow rate 	[Sm^3/d]
%	u(2*(k-1)+2): W(k).T   	Inlet W temperature                	[K]
%	u(2*nW+1):    P        	Mixer node pressure                	[kPa]
%
%   Outputs:
%   y(1): Qw        		Water brine std. vol. flow rate    	[Sm^3/d]
%   y(2): T         		Temperature                        	[K]
%
%   Initial conditions set:
%   There are no initial conditions for this model.

%   2023. Coded by Diego Trica (diegotrica@gmail.com)

%Blackoil object class property package. On the main script the
%user must declare PP as global and PP=blackoil.
global PP

%Setting model parameters:
Cws = par{1};				%Water salt concentration			[g/L]
nW 	= par{2};   			%Number of W inlet streams 			[]

%Water salt mass fraction []
ws=Cws/(PP.rho_H2Ostd+Cws);

%Water brine mass density [kg/m^3]
rWstd=WaterStdMassDens(PP,ws);

%Water specific heat capacity [kJ/kg K]
cPw=WaterSpecHeatCap(PP,PP.Pref,PP.Tref,ws);

switch flag,
%----- S-FUNCTION BUILDING PARAMETERS -------------------------------
    case 0,
sys(1)=0;%Number of continuous states
sys(2)=0;%Number of discrete states
sys(3)=2;%Number of outputs
sys(4)=2*nW+1;%Number of inputs
sys(5)=0;%Reserved for root finding. Must be zero.
sys(6)=1;%Direct Feedthrough S-function parameter
sys(7)=1;%Number of sample times (at least one sample time is needed)

%----- States initial condition:
x0=[];

%----- Other building parameters:
str=[];     %str is always an empty matrix
ts=[0 0];   %Sample time array [0 0] for continuous system
simStateCompliance='UnknownSimState';
%--------------------------------------------------------------------
%----- STATE DERIVATIVES --------------------------------------------
    case 1,
sys=[];
%--------------------------------------------------------------------
    case 2,
sys=[];
%----- OUPUTS -------------------------------------------------------
    case 3,
%----- Inputs:
for k=1:nW
W(k).Qw		= u(2*(k-1)+1);
W(k).T     	= u(2*(k-1)+2);
end
P       	= u(2*nW+1);

%----- Output equations:
%Water brine std. vol. balance [Sm^3/d]
Qw = sum([W.Qw]);

%Outlet temperature [K]
T = PP.Tref + ...
	sum([W.Qw]*rWstd*cPw.*([W.T]-PP.Tref))/(Qw*rWstd*cPw);

%----- Outputs for sys output flag=3:
sys(1)=Qw;
sys(2)=T;
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