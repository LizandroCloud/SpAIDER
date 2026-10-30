function [sys,x0,str,ts,simStateCompliance]=mixer_comp(t,x,u,flag,css,par)
%mixer3_gas - Gas stream inlet mixer for unit compression non-linear
%   model S-Function
%
%   This model describes the dynamic behavior of a mixer with three 
%	gas-only inlet streams for usage in unit compression system of
%   oil & gas upstream process, using the black-oil approach.
%
%   [sys,x0,str,ts,simStateCompliance]=mixer_comp(t,x,u,flag,css,par)
%
%   sys is the S-function output, which depends on which flag is 
%   being inputted by simulink. x0 are the states initial condi-
%   tions. See the S-function Level-1 MATLAB reference page for more
%   details about str, ts and simStateCompliance S-Function ouputs.
%
%   t is the independent variable time. x are the state variables.
%   u are the input variables, combined into one vector of manipula-
%   ted and disturbs. flag is the Simulink variable inputted to get a
%   specific output from the S-function. css are the fixed set passed
%   by the user to compute x0 for steady-state condition. par are the
%   model parameters.
%
%   Model parameters:
%   par{1}: d60     Oil specific gravity                     []
%   par{2}: dg      Gas specific gravity                     []
%
%   States:
%   There are no states for this model.
%
%   Inputs:
%	u(1): Qg_sg     Rec. surge stream std. vol. flow rate    [Sm^3/d]
%   u(2): T_sg      Rec. surge stream temperature            [K]
%	u(3): Qg_comp   Prev. comp. stage std. vol. flow rate    [Sm^3/d]
%   u(4): T_comp    Prev. comp. stage temperature            [K]
%   u(5): T_feed    Feed gas stream temperature              [K]
%   u(6): Qg_feed   Feed stream gas std. vol. flow rate      [Sm^3/d]
%   u(7): P         Mixer node pressure                      [kPa]
%
%   Outputs:
%   y(1): Qg_out    Outlet gas std. vol. flow rate  		 [Sm^3/d]
%   y(2): T         Outlet temperature                       [K]
%
%   Fixed set for steady-state condition:
%	The model has no dynamics, for which the call with flag=3
%	computes a steady-state condition.

%   2022. Coded by Diego Trica (diegotrica@gmail.com)

%Blackoil object class property package. On the main script the
%user must declare PP as global and PP=blackoil.
    global PP

%Setting model parameters:
    d60         = par{1};   %Oil specific gravity []    
    dg          = par{2};   %Gas specific gravity []
   
%Oil specific heat capacity [kJ/kg K]
%cPo=OilSpecHeatCap(PP,d60,Tbn);

%Gas specific heat capacity [kJ/kg K]
cPg=GasSpecHeatCap(PP,dg);

%Gas liberation specific enthapy [kJ/kg]
DeOG=GasLatHeat(PP,d60,dg);

switch flag,
%----- S-FUNCTION BUILDING PARAMETERS -------------------------------
    case 0,
sys(1)=0;%Number of continuous states
sys(2)=0;%Number of discrete states
sys(3)=2;%Number of outputs
sys(4)=7;%Number of inputs
sys(5)=0;%Reserved for root finding. Must be zero.
sys(6)=1;%Direct Feedthrough S-function parameter
sys(7)=1;%Number of sample times (at least one sample time is needed)

%----- States at steady-state condition:
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
Qg_sg	= u(1);
T_sg    = u(2);
Qg_comp	= u(3);
T_comp  = u(4);
Qg_feed = u(5);
T_feed  = u(6);
P       = u(7);

%----- Model internal variables:
%Inlet sg, comp, feed gas residual specific enthalpy [kJ/kg K]
eGres_sg=GasSpecResEnth(PP,P,T_sg,dg);
eGres_comp=GasSpecResEnth(PP,P,T_comp,dg);
eGres_feed=GasSpecResEnth(PP,P,T_feed,dg);

%----- Output equations:
%Gas std. vol. balance
%Qg_feed=Qg_out-Qg_sg-Qg_comp;
Qg_out=Qg_feed+Qg_sg+Qg_comp;

%Energy balance: isenthalpic expansion - Iterative calculation
T=T_feed; iter=0; maxiter=100; rtol=1e-5;
while iter<maxiter
    iter=iter+1;
    T0=T;
	
%Outlet gas specific residual enthalpy [kJ/kg]
eGres=GasSpecResEnth(PP,P,T,dg);

%Gas specific residual enthalpy derivative do T [kJ/kg K]
deGresdT=dGasSpecResEnthdT(PP,P,T,dg);

%Outlet temperature [K]
f=...
 +Qg_out*PP.rho_AIRstd*dg*(eGres+DeOG+cPg*(T-PP.Tref)) ...
 -Qg_sg*PP.rho_AIRstd*dg*(eGres_sg+DeOG+cPg*(T_sg-PP.Tref)) ...
 -Qg_comp*PP.rho_AIRstd*dg*(eGres_comp+DeOG+cPg*(T_comp-PP.Tref)) ...
 -Qg_feed*PP.rho_AIRstd*dg*(eGres_feed+DeOG+cPg*(T_feed-PP.Tref));
dfdT=Qg_out*PP.rho_AIRstd*dg*(deGresdT+cPg);
T=T-f/dfdT;
	rdev=abs((T-T0)/T0);
	if rdev<rtol
		break
	end
end

%----- Outputs for sys output flag=3:
sys(1)=Qg_out;
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