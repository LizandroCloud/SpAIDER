function [sys,x0,str,ts,simStateCompliance]=hexch_liq(t,x,u,flag,cinit,par)
%hexch_liq - Liquid stream heat exchanger non-linear model S-Function
%
%   This model describes the dynamic behavior of a heat exchanger on
%   counter-current arrangement with two inlet triphase oil-gas-water
%   phase streams for oil & gas upstream process, using the black-oil
%   approach.
%
%   [sys,x0,str,ts,simStateCompliance]=hexch_liq(t,x,u,flag,cinit,par)
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
%   par{5}: UA      Global heat transfer coefficient         	[kW/K]
%
%   States:
%   There are no states for this model.
%
%   Inputs:
%	u(1):  Qf_inC   Inlet cold side liq. std. vol. flow rate [Sm^3/d]
%	u(2):  Qg_inC   Inlet cold side gas std. vol. flow rate  [Sm^3/d]
%	u(3):  BSW_inC  Inlet cold side basic sediments & water  []
%   u(4):  P_inC    Inlet cold side pressure                 [kPa]
%   u(5):  T_inC    Inlet cold side temperature              [K]
%   u(6):  P_outC   Outlet cold side pressure                [kPa]
%	u(7):  Qf_inH   Inlet hot side liq. std. vol. flow rate  [Sm^3/d]
%	u(8):  Qg_inH   Inlet hot side gas std. vol. flow rate   [Sm^3/d]
%	u(9):  BSW_inH  Inlet hot side basic sediments & water   []
%   u(10): P_inH    Inlet hot side pressure                  [kPa]
%   u(11): T_inH    Inlet hot side temperature               [K]
%   u(12): P_outH   Outlet hot side pressure                 [kPa]
%
%   Outputs:
%   y(1): Qf_C      Out. cold side liq. std. vol. flow rate  [Sm^3/d]
%   y(2): Qg_C      Outlet cold side gas std. vol. flow rate [Sm^3/d]
%   y(3): BSW_C     Outlet cold side basic sediments & water []
%   y(4): T_C       Outlet cold side temperature             [K]
%   y(5): Qf_H      Out. hot side liq. std. vol. flow rate   [Sm^3/d]
%   y(6): Qg_H      Outlet hot side gas std. vol. flow rate  [Sm^3/d]
%   y(7): BSW_H     Outlet hot side basic sediments & water  []
%   y(8): T_H       Outlet hot side temperature              [K]
%
%   Initial conditions set:
%   There are no initial conditions for this model.

%   2022. Coded by Diego Trica (diegotrica@gmail.com)

%Blackoil object class property package. On the main script the
%user must declare PP as global and PP=blackoil.
global PP

%Setting model parameters:
h           = 1e-2;     %regularization function parameter  []
d60         = par{1};   %Oil specific gravity []
Tbn         = par{2};   %Oil normal boiling point [K]
Cws     	= par{3};	%Water salt concentration [g/L]
dg          = par{4};   %Gas specific gravity []
UA          = par{5};   %Global heat transfer coefficient [kW/K]
  
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

switch flag,
%----- S-FUNCTION BUILDING PARAMETERS -------------------------------
    case 0,
sys(1)=0;%Number of continuous states
sys(2)=0;%Number of discrete states
sys(3)=8;%Number of outputs
sys(4)=12;%Number of inputs
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
nIn=6;

Qf_inC	= u(1);
Qg_inC	= u(2);
BSW_inC	= u(3);
P_inC   = u(4);
T_inC   = u(5);
P_outC  = u(6);

Qf_inH	= u(nIn+1);
Qg_inH	= u(nIn+2);
BSW_inH	= u(nIn+3);
P_inH   = u(nIn+4);
T_inH   = u(nIn+5);
P_outH  = u(nIn+6);

%----- Model internal variables:
%Cold/Hot streams gas spec. res. enthalpy derivative to T [kJ/kg K]
deGresdT_C=dGasSpecResEnthdT(PP,P_inC,T_inC,dg);
deGresdT_H=dGasSpecResEnthdT(PP,P_inH,T_inH,dg);

%Cold/Hot streams heat capacity rate [kW/K]
Ccold=(Qf_inC*((1-BSW_inC)*d60*PP.rho_H2Ostd*cPo+BSW_inC*rWstd*cPw)...
 +Qg_inC*dg*PP.rho_AIRstd*(deGresdT_C+cPg))/(24*3600);
Chot=(Qf_inH*((1-BSW_inH)*d60*PP.rho_H2Ostd*cPo+BSW_inH*rWstd*cPw)...
 +Qg_inH*dg*PP.rho_AIRstd*(deGresdT_H+cPg))/(24*3600);

%Heat capacity rate [kW/K]
Cmin=0.5*(1+tanh((Ccold-Chot)/h))*Chot+0.5*(1+tanh((Chot-Ccold)/h))*Ccold;
Cmax=0.5*(1+tanh((Ccold-Chot)/h))*Ccold+0.5*(1+tanh((Chot-Ccold)/h))*Chot;

%Thermodynamic limit for heat exchange
Qmax=Cmin*(T_inH-T_inC);

%Number of transfer units
NTU=UA/Cmin;

%Heat capacity rate ratio
Cr=Cmin/Cmax;
 
%Heat transfer effectiveness: Counter-current flow 
eff=(1-exp(-NTU*(1-Cr)))/(1-Cr*exp(-NTU*(1-Cr)));

%Heat transferred effectiveness
Qexch=eff*Qmax;

%----- Output equations:
uCold=[Qf_inC; Qg_inC; BSW_inC; P_inC; T_inC; P_outC; Qexch];
yCold=heater_liq([],[],uCold,3,[],{par{1},par{2},par{3},par{4}});

uHot=[Qf_inH; Qg_inH; BSW_inH; P_inH; T_inH; P_outH; Qexch];
yHot=cooler_liq([],[],uHot,3,[],{par{1},par{2},par{3},par{4}});
%----- Outputs for sys output flag=3:
sys=[yCold;yHot];
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