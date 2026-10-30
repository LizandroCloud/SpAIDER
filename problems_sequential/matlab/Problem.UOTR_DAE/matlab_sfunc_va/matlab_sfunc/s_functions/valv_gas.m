function [sys,x0,str,ts,simStateCompliance]=valv_gas(t,x,u,flag,cinit,par)
%valv_gas - Gas flow valve non-linear model S-Function
%
%   This model describes the dynamic behavior of a valve with a gas-
%   only inlet/outlet streams for oil & gas upstream processes, using
%   the black-oil approach.
%
%   [sys,x0,str,ts,simStateCompliance]=valv_gas(t,x,u,flag,cinit,par)
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
%   par{2}: Cv      Valve coefficient constant US units      [gpm]
%   par{3}: Fgm     Specific heat capacity ratio factor      []
%   par{4}: DPrT    Pressure drop ratio factor               []
%
%   States:
%   There are no states for this model.
%
%   Inputs:
%   u(1): P_in      Inlet pressure                           [kPa]
%   u(2): T_in      Inlet temperature                        [K]
%   u(3): P_out     Outlet pressure                          [kPa]
%   u(4): c         Valve stem position                      []
%
%   Outputs:
%   y(1): Qg        Gas std. vol. flow rate                  [Sm^3/d]
%   y(2): T         Temperature                              [K]
%
%   Initial conditions set:
%   There are no initial conditions for this model.

%   2022. Coded by Diego Trica (diegotrica@gmail.com)

%Blackoil object class property package. On the main script the
%user must declare PP as global and PP=blackoil.
global PP

%Setting model parameters:
dg          = par{1};   %Gas specific gravity []
Cv          = par{2};   %Valve coeff. constant US units [gpm]
Fgm         = par{3};   %Specific heat capacity ratio factor []
DPrT        = par{4};   %Pressure drop ratio factor []

%Gas specific heat capacity [kJ/kg K]
cPg=GasSpecHeatCap(PP,dg);
   
%Valve coefficient constant [m^2]
Av=Cv*(0.0000630902)*sqrt(PP.rho_H2Ostd*14.696/101325);

%Pressure drop ratio at choked flow []
DPrchok = Fgm*DPrT;

switch flag,
%----- S-FUNCTION BUILDING PARAMETERS -------------------------------
    case 0,
sys(1)=0;   %Number of continuous states
sys(2)=0;   %Number of discrete states
sys(3)=2;   %Number of outputs
sys(4)=4;   %Number of inputs
sys(5)=0;   %Reserved for root finding. Must be zero.
sys(6)=1;   %Direct Feedthrough S-function parameter
sys(7)=1;   %Number of sample times (only for discrete-time stystem)

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
P_in    = u(1);
T_in	= u(2);
P_out   = u(3);
c       = u(4);

%----- Model internal variables:
%Inlet gas specific residual enthalpy [kJ/kg]
eGres_in=GasSpecResEnth(PP,P_in,T_in,dg);

%Inlet gas compressibility factor []
Zg_in=GasCompFactor(PP,P_in,T_in,dg);

%Inlet gas mass density [kg/m^3]
rG_in=P_in*(dg*PP.MM_Air)/(PP.R*Zg_in*T_in);

%Pressure drop [kPa]
DP=P_in-P_out;

%Pressure drop ratio []
DPr=DP/P_in;

%Linear inherent characteristic function
fc=c;

%----- Output equations:
%Pressure-flow relation
DPrflow=0.5*(1+tanh(((DPrchok-DPr)+1e-3)/1e-4))*DPr...
 + 0.5*(1+tanh((DPr-DPrchok)/1e-4))*DPrchok;
Fexp=1-DPrflow/(3*DPrchok); %Expansion factor []
Qg=(24*3600)/(dg*PP.rho_AIRstd)*fc*Av*Fexp*sqrt(DPrflow*P_in*rG_in*1e3);

%Energy balance: isenthalpic expansion - Iterative calculation
T=T_in; iter=0; maxiter=100; rtol=1e-5;
while iter<maxiter
    iter=iter+1;
    T0=T;
    eGres=GasSpecResEnth(PP,P_out,T,dg);    %Outlet gas specific residual enthalpy [kJ/kg]
    T=T_in + (eGres_in-eGres)/cPg;
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