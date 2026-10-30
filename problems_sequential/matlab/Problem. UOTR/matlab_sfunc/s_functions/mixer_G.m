function [sys,x0,str,ts,simStateCompliance]=mixer_G(t,x,u,flag,cinit,par)
%mixer_gas - Gas streams mixer non-linear model S-Function
%
%   This model describes the dynamic behavior of a mixer with gas
%   only inlet streams for oil & gas upstream process, using the 
%	black-oil approach.
%
%   [sys,x0,str,ts,simStateCompliance]=mixer_G(t,x,u,flag,cinit,par)
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
%   par{1}: dg      Gas specific gravity                     []
%   par{2}: nG      Number of G inlet streams                []
%
%   States:
%   There are no states for this model.
%
%   Inputs:
%	u(2*(k-1)+1): G(k).Qg   Inlet G gas std. vol. flow rate 	[Sm^3/d]
%   u(2*(k-1)+2): G(k).T    Inlet G temperature                 [K]
%	u(2*nG+1):    P        	Mixer node pressure                	[kPa]
%
%   Outputs:
%   y(1): Qg        		Gas std. vol. flow rate    			[Sm^3/d]
%   y(2): T         		Temperature                        	[K]
%
%   Initial conditions set:
%   There are no initial conditions for this model.

%   2023. Coded by Diego Trica (diegotrica@gmail.com)

%Blackoil object class property package. On the main script the
%user must declare PP as global and PP=blackoil.
global PP

%Setting model parameters:
dg  = par{1};   %Gas specific gravity 				[]
nG 	= par{2};   %Number of G inlet streams 			[]

%Gas specific heat capacity [kJ/kg K]
cPg=GasSpecHeatCap(PP,dg);

switch flag,
%----- S-FUNCTION BUILDING PARAMETERS -------------------------------
    case 0,
sys(1)=0;%Number of continuous states
sys(2)=0;%Number of discrete states
sys(3)=2;%Number of outputs
sys(4)=2*nG+1;%Number of inputs
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
for k=1:nG
G(k).Qg		= u(2*(k-1)+1);
G(k).T     	= u(2*(k-1)+2);
end
P       	= u(2*nG+1);

%----- Model internal variables:
for k=1:nG
%Inlet 1,2 gas residual specific enthalpy [kJ/kg K]
G(k).eGres=GasSpecResEnth(PP,P,G(k).T,dg);
end

%----- Output equations:
%Gas std. vol. balance
Qg=sum([G.Qg]);

%Energy balance: isenthalpic expansion - Iterative calculation
T=min([G.T]); iter=0; maxiter=100; rtol=1e-5;
while iter<maxiter
    iter=iter+1;
    T0=T;
	
%Outlet gas specific residual enthalpy [kJ/kg]
eGres=GasSpecResEnth(PP,P,T,dg);

%Gas specific residual enthalpy derivative do T [kJ/kg K]
deGresdT=dGasSpecResEnthdT(PP,P,T,dg);

%Outlet temperature [K]
f = sum([G.Qg]*PP.rho_AIRstd*dg.*([G.eGres]+cPg.*([G.T]-PP.Tref))) ...
	- Qg*PP.rho_AIRstd*dg*(eGres+cPg*(T-PP.Tref));
dfdT = - Qg*PP.rho_AIRstd*dg*(deGresdT+cPg);
T=T-f/dfdT;

%{
T = PP.Tref + (...
	+ sum([G.Qg]*PP.rho_AIRstd*dg.*([G.eGres]+cPg.*([G.T]-PP.Tref))) ...
	- Qg*PP.rho_AIRstd*dg*(eGres))/(Qg*PP.rho_AIRstd*dg*cPg);
%}
	rdev=abs((T-T0)/T0);
	if rdev<rtol
		break
	end
end

%----- Outputs for sys output flag=3:
sys(1)=Qg;
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