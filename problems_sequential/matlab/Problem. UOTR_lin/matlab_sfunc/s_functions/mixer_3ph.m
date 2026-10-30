function [sys,x0,str,ts,simStateCompliance]=mixer_3ph(t,x,u,flag,css,par)
%mixer_3ph - 3-phase streams mixer non-linear model S-Function
%
%   This model describes the dynamic behavior of mixer with two inlet
%   triphase oil-gas-water stream for oil & gas upstream process, 
%   using the black-oil approach.
%
%   [sys,x0,str,ts,simStateCompliance]=mixer_3ph(t,x,u,flag,css,par)
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
%   par{2}: Tbn     Oil normal boiling point                 [K]
%   par{3}: dg      Gas specific gravity                     []
%
%   States:
%   There are no states for this model.
%
%   Inputs:
%	u(1): Qf_1      Inlet 1 liquid std. vol. flow rate       [Sm^3/d]
%	u(2): Qg_1      Inlet 1 gas std. vol. flow rate          [Sm^3/d]
%	u(3): BSW_1     Inlet 1 basic sediments & water          []
%	u(4): ws_1      Inlet 1 water salt mass fraction         []
%   u(5): T_1       Inlet 1 temperature                      [K]
%	u(6): Qf_2      Inlet 2 liquid std. vol. flow rate       [Sm^3/d]
%	u(7): Qg_2      Inlet 2 gas std. vol. flow rate          [Sm^3/d]
%	u(8): BSW_2     Inlet 2 basic sediments & water          []
%	u(9): ws_2      Inlet 2 water salt mass fraction         []
%   u(10): T_2      Inlet 2 temperature                      [K]
%   u(11): P        Mixer node pressure                      [kPa]
%
%   Outputs:
%   y(1): Qf        Liquid std. vol. flow rate               [Sm^3/d]
%   y(2): Qg        Gas std. vol. flow rate                  [Sm^3/d]
%   y(3): BSW       Basic sediments & water                  []
%   y(4): ws        Water salt mass fraction                 []
%   y(5): T         Temperature                              [K]
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
	Tbn         = par{2};   %Oil normal boiling point [K]
    dg          = par{3};   %Gas specific gravity []
   
%Oil specific heat capacity [kJ/kg K]
cPo=OilSpecHeatCap(PP,d60,Tbn);

%Gas specific heat capacity [kJ/kg K]
cPg=GasSpecHeatCap(PP,dg);

%Gas liberation specific enthapy [kJ/kg]
DeOG=GasLatHeat(PP,d60,dg);

switch flag,
%----- S-FUNCTION BUILDING PARAMETERS -------------------------------
    case 0,
sys(1)=0;%Number of continuous states
sys(2)=0;%Number of discrete states
sys(3)=5;%Number of outputs
sys(4)=11;%Number of inputs
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
Qf_1	= u(1);
Qg_1	= u(2);
BSW_1	= u(3);
ws_1	= u(4);
T_1     = u(5);
Qf_2	= u(6);
Qg_2	= u(7);
BSW_2	= u(8);
ws_2	= u(9);
T_2     = u(10);
P       = u(11);

%----- Model internal variables:
%Inlet 1,2 solution gas/oil ratio [Sm^3/Sm^3]
Rs_1=SolGasOilRatio(PP,P,T_1,d60,dg);
Rs_2=SolGasOilRatio(PP,P,T_2,d60,dg);

%Inlet 1,2 gas specific residual enthalpy [kJ/kg]
eGres_1=GasSpecResEnth(PP,P,T_1,dg);
eGres_2=GasSpecResEnth(PP,P,T_1,dg);

%Inlet 1,2 water brine mass density std. cond. [kg/m^3]
rWstd_1=WaterStdMassDens(PP,ws_1);
rWstd_2=WaterStdMassDens(PP,ws_2);

%Inlet 1,2 water brine specific heat capacity [kJ/kg K]
cPw_1=WaterSpecHeatCap(PP,PP.Pref,PP.Tref,ws_1);
cPw_2=WaterSpecHeatCap(PP,PP.Pref,PP.Tref,ws_2);

%----- Output equations:
%Liquid std. vol. balance
Qf=Qf_1+Qf_2;

%Basic sediments & water
BSW=(Qf_1*BSW_1+Qf_2*BSW_2)/Qf;

%Water salt mass fraction [] - Iterative calculation
ws=(Qf_1*BSW_1*rWstd_1*ws_1+Qf_2*BSW_2*rWstd_2*ws_2)/(Qf*BSW*1030);
iter=0; maxiter=100; rtol=1e-5;
while iter<maxiter
    iter=iter+1;
    ws0=ws;
    rWstd=WaterStdMassDens(PP,ws); %Water brine mass density [kg/m^3]    
    ws=(Qf_1*BSW_1*rWstd_1*ws_1+Qf_2*BSW_2*rWstd_2*ws_2)/(Qf*BSW*rWstd);
    rdev=abs((ws-ws0)/ws0);
    if rdev<rtol
        break
    end
end

%Water brine specific heat capacity [kJ/kg K]
cPw=WaterSpecHeatCap(PP,PP.Pref,PP.Tref,ws);

%Energy balance: isenthalpic expansion - Iterative calculation
T=T_1; iter=0; maxiter=100; rtol=1e-5;
while iter<maxiter
    iter=iter+1;
    T0=T;
    
%Solution gas/oil ratio [Sm^3/Sm^3]
Rseq=SolGasOilRatio(PP,P,T,d60,dg);

%Sum of inlet gas std. vol.
Qg_sumIn=Qf_1*(1-BSW_1)*Rs_1+Qf_2*(1-BSW_2)*Rs_2+Qg_1+Qg_2;

%Gas std. vol. unbalance
Qg_unbal=Qf*(1-BSW)*Rseq-Qg_sumIn;

%Std. vol. balance for liberated gas [Sm^3/d]
Qg=0.5*(1+tanh(-(Qg_unbal)/1e-6))*(-Qg_unbal);

%Outlet gas specific residual enthalpy [kJ/kg]
eGres=GasSpecResEnth(PP,P,T,dg);

%Gas specific residual enthalpy derivative do T [kJ/kg K]
deGresdT=dGasSpecResEnthdT(PP,P,T,dg);

%Outlet temperature [K]
f=Qf*((1-BSW)*d60*PP.rho_H2Ostd*cPo*(T-PP.Tref)+BSW*rWstd*cPw*(T-PP.Tref))...
 + Qg*PP.rho_AIRstd*dg*(eGres+DeOG+cPg*(T-PP.Tref)) ...
 - Qf_1*((1-BSW_1)*d60*PP.rho_H2Ostd*cPo*(T_1-PP.Tref)+BSW_1*rWstd_1*cPw_1*(T_1-PP.Tref))...
 - Qg_1*PP.rho_AIRstd*dg*(eGres_1+DeOG+cPg*(T_1-PP.Tref))...
 - Qf_2*((1-BSW_2)*d60*PP.rho_H2Ostd*cPo*(T_2-PP.Tref)+BSW_2*rWstd_2*cPw_2*(T_2-PP.Tref))...
 - Qg_2*PP.rho_AIRstd*dg*(eGres_2+DeOG+cPg*(T_2-PP.Tref));
dfdT=Qf*((1-BSW)*d60*PP.rho_H2Ostd*cPo+BSW*rWstd*cPw)...
 + Qg*PP.rho_AIRstd*dg*(deGresdT+cPg);
T=T-f/dfdT;
	rdev=abs((T-T0)/T0);
	if rdev<rtol
		break
	end
end

%----- Outputs for sys output flag=3:
sys(1)=Qf;
sys(2)=Qg;
sys(3)=BSW;
sys(4)=ws;
sys(5)=T;
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