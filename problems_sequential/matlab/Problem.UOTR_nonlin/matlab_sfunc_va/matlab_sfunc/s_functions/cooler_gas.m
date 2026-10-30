function [sys,x0,str,ts,simStateCompliance]=cooler_gas(t,x,u,flag,cinit,par)
%cooler_gas - Gas stream heater non-linear model S-Function
%
%   This model describes the dynamic behavior of a cooler with a gas-
%	only inlet stream for oil & gas upstream process, using the
%	black-oil approach.
%
%   [sys,x0,str,ts,simStateCompliance]=cooler_gas(t,x,u,flag,cinit,par)
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
%
%   States:
%   There are no states for this model.
%
%   Inputs:
%   u(1): Qg_in     Inlet std. vol. flow rate                [Sm^3/d]
%	u(2): P_in      Inlet pressure                           [kPa]
%   u(3): T_in      Inlet temperature                        [K]
%   u(4): P_out     Outlet pressure                          [kPa]
%   u(5): Q_out     Heat removed from the system             [kW]
%
%   Outputs:
%   y(1): Qg        Gas std. vol. flow rate                  [Sm^3/d]
%   y(2): T         Temperature                              [K]
%
%   Initial conditions set:
%   There are no initial conditions for this model.

%   2023. Coded by Diego Trica (diegotrica@gmail.com)

%Blackoil object class property package. On the main script the
%user must declare PP as global and PP=blackoil.
global PP

%Setting model parameters:
dg = par{1};   %Gas specific gravity []
   
%Gas specific heat capacity [kJ/kg K]
cPg=GasSpecHeatCap(PP,dg);

%Gas liberation specific enthapy [kJ/kg]
%DeOG=GasLatHeat(PP,d60,dg);

switch flag,
%----- S-FUNCTION BUILDING PARAMETERS -------------------------------
    case 0,
sys(1)=0;%Number of continuous states
sys(2)=0;%Number of discrete states
sys(3)=2;%Number of outputs
sys(4)=5;%Number of inputs
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
Qg_in   = u(1);
P_in    = u(2);
T_in	= u(3);
P_out   = u(4);
Q_out   = u(5);

%----- Model internal variables:
%Inlet gas specific residual enthalpy [kJ/kg]
eGres_in=GasSpecResEnth(PP,P_in,T_in,dg);

%----- Output equations:
%Gas std. vol. balance [Sm^3/d]
Qg=Qg_in;

%Energy balance: isenthalpic expansion - Iterative calculation
T=T_in; iter=0; maxiter=100; rtol=1e-5;
while iter<maxiter
    iter=iter+1;
    T0=T;

%Outlet gas specific residual enthalpy [kJ/kg]
eGres = GasSpecResEnth(PP,P_out,T,dg);

%Outlet temperature [K]
%{
T = T_in + (eGres_in - eGres)/cPg ...
	- Q_out*(24*3600)/(Qg*PP.rho_AIRstd*dg*cPg);
%}
f = Qg/(24*3600)*PP.rho_AIRstd*dg*(eGres_in+cPg*(T_in-PP.Tref)) ...
	- Qg/(24*3600)*PP.rho_AIRstd*dg*(eGres+cPg*(T-PP.Tref)) ...
	- Q_out;
deGresdT = dGasSpecResEnthdT(PP,P_out,T,dg);
dfdT = -Qg/(24*3600)*PP.rho_AIRstd*dg*(deGresdT+cPg);
T = T - f/dfdT;

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