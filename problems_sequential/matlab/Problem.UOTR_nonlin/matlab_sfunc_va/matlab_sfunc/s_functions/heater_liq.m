function [sys,x0,str,ts,simStateCompliance]=heater_liq(t,x,u,flag,cinit,par)
%heater_liq - Liquid stream heater non-linear model S-Function
%
%   This model describes the dynamic behavior of a heater with a tri-
%	phase oil-gas-water inlet stream for oil & gas upstream process, 
%   using the black-oil approach.
%
%   [sys,x0,str,ts,simStateCompliance]=heater_liq(t,x,u,flag,cinit,par)
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
%   par{1}: d60     Oil specific gravity                     	[]
%   par{2}: Tbn     Oil normal boiling point                 	[K]
%   par{3}: Cws     Water salt concentration                	[g/L]
%   par{4}: dg      Gas specific gravity                     	[]
%
%   States:
%   There are no states for this model.
%
%   Inputs:
%	u(1): Qf_in     Inlet liquid std. vol. flow rate         [Sm^3/d]
%	u(2): Qg_in     Inlet gas std. vol. flow rate            [Sm^3/d]
%	u(3): BSW_in    Inlet basic sediments & water            []
%   u(4): P_in      Inlet pressure                           [kPa]
%   u(5): T_in      Inlet temperature                        [K]
%   u(6): P_out     Outlet pressure                          [kPa]
%   u(7): Q_in      Heat applied to the system               [kW]
%
%   Outputs:
%   y(1): Qf        Liquid std. vol. flow rate               [Sm^3/d]
%   y(2): Qg        Gas std. vol. flow rate                  [Sm^3/d]
%   y(3): BSW       Basic sediments & water                  []
%   y(4): T         Temperature                              [K]
%
%   Initial conditions set:
%   There are no initial conditions for this model.

%   2022. Coded by Diego Trica (diegotrica@gmail.com)

%Blackoil object class property package. On the main script the
%user must declare PP as global and PP=blackoil.
global PP

%Setting model parameters:
d60         = par{1};   %Oil specific gravity 				[]
Tbn         = par{2};   %Oil normal boiling point 			[K]
Cws     	= par{3};	%Water salt concentration			[g/L]	
dg          = par{4};   %Gas specific gravity 				[]

%Water salt mass fraction []
ws=Cws/(999.054+Cws);
   
%Water brine mass density [kg/m^3]
rWstd=WaterStdMassDens(PP,ws);   
   
%Oil specific heat capacity [kJ/kg K]
cPo=OilSpecHeatCap(PP,d60,Tbn);

%Water specific heat capacity [kJ/kg K]
cPw=WaterSpecHeatCap(PP,PP.Pref,PP.Tref,ws);

%Gas specific heat capacity [kJ/kg K]
cPg=GasSpecHeatCap(PP,dg);

%Gas liberation specific enthapy [kJ/kg]
DeOG=GasLatHeat(PP,d60,dg);

switch flag,
%----- S-FUNCTION BUILDING PARAMETERS -------------------------------
    case 0,
sys(1)=0;%Number of continuous states
sys(2)=0;%Number of discrete states
sys(3)=4;%Number of outputs
sys(4)=7;%Number of inputs
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
%----- State derivatives for sys output flag=1
sys=[];
%--------------------------------------------------------------------
    case 2,
sys=[];
%----- OUPUTS -------------------------------------------------------
    case 3,
%----- Inputs:
Qf_in	= u(1);
Qg_in	= u(2);
BSW_in	= u(3);
P_in    = u(4);
T_in	= u(5);
P_out   = u(6);
Q_in    = u(7);

%----- Model internal variables:
%Inlet solution gas/oil ratio [Sm^3/Sm^3]
Rs_in=SolGasOilRatio(PP,P_in,T_in,d60,dg);

%Inlet gas specific residual enthalpy [kJ/kg]
eGres_in=GasSpecResEnth(PP,P_in,T_in,dg);

%----- Output equations:
%Liquid std. vol. balance
Qf=Qf_in;

%Water std. vol. balance
BSW=BSW_in;

%Energy balance: isenthalpic expansion - Iterative calculation
T=T_in; iter=0; maxiter=100; rtol=1e-5;
while iter<maxiter
    iter=iter+1;
    T0=T;
    
%Solution gas/oil ratio [Sm^3/Sm^3]
Rseq = SolGasOilRatio(PP,P_out,T,d60,dg);

%Std. vol. balance for liberated gas [Sm^3/d]
Qg=Qg_in+0.5*(1+tanh(-(Rseq-Rs_in)/1e-6))*Qf*(1-BSW)*(Rs_in-Rseq);

%Outlet gas specific residual enthalpy [kJ/kg]
eGres=GasSpecResEnth(PP,P_out,T,dg);

%Outlet temperature [K]
%{
T = PP.Tref + (...
	Q_in*(24*3600)...
	+ Qf*(1-BSW)*d60*PP.rho_H2Ostd*cPo*(T_in-PP.Tref) ...
	+ Qf*BSW*rWstd*cPw*(T_in-PP.Tref) ...
	- Qg*PP.rho_AIRstd*dg*(eGres+DeOG))/(...
		Qf*((1-BSW)*d60*PP.rho_H2Ostd*cPo + BSW*rWstd*cPw) + ...
		Qg*PP.rho_AIRstd*dg*cPg);
%}
T=T_in...
 +Q_in/(Qf/(24*3600)*((1-BSW)*d60*PP.rho_H2Ostd*cPo+BSW*rWstd*cPw))...
 +Qg_in*dg*PP.rho_AIRstd/(Qf*((1-BSW)*d60*PP.rho_H2Ostd*cPo+BSW*rWstd*cPw))*(eGres_in+DeOG+cPg*(T_in-PP.Tref))...
 -Qg*dg*PP.rho_AIRstd/(Qf*((1-BSW)*d60*PP.rho_H2Ostd*cPo+BSW*rWstd*cPw))*(eGres+DeOG+cPg*(T-PP.Tref));

	rdev=abs((T-T0)/T0);
	if rdev<rtol
		break
	end
end

%----- Outputs for sys output flag=3:
sys(1)=Qf;
sys(2)=Qg;
sys(3)=BSW;
sys(4)=T;
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