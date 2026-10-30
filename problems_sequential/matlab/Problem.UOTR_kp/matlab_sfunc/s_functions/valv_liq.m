function [sys,x0,str,ts,simStateCompliance]=valv_liq(t,x,u,flag,cinit,par)
%valv_liq - Liquid flow valve non-linear model S-Function
%
%   This model describes the dynamic behavior of a valve with a li-
%   quid oil-water inlet stream for oil & gas upstream processes, 
%   using the black-oil approach.
%
%   [sys,x0,str,ts,simStateCompliance]=valv_liq(t,x,u,flag,cinit,par)
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
%   par{1}: d60     Oil specific gravity                     []
%   par{2}: Tbn     Oil normal boiling point                 [K]
%   par{3}: Cws     Water salt concentration                 [g/L]
%   par{4}: dg      Gas specific gravity                     []
%   par{5}: Cv      Valve constant US units      			 [gpm]
%   par{6}: FL      Liquid pressure recovery factor          []
%
%   States:
%   There are no states for this model.
%
%   Inputs:
%   u(1): BSW_in    Inlet basic sediments & water            []
%   u(2): P_in      Inlet pressure                           [kPa]
%   u(3): T_in      Inlet temperature                        [K]
%   u(4): P_out     Outlet pressure                          [kPa]
%   u(5): c     	Valve stem position             		 []
%
%   Outputs:
%   y(1): Qf        Liquid std. vol. flow rate               [Sm^3/d]
%   y(2): Qg        Gas std. vol. flow rate                  [Sm^3/d]
%   y(3): BSW       Basic sediments & water                  []
%   y(4): T         Temperature                              [K]
%
%   Initial conditions set:
%   There are no initial conditions for this model.

%   2023. Coded by Diego Trica (diegotrica@gmail.com)

%Blackoil object class property package. On the main script the
%user must declare PP as global and PP=blackoil.
global PP

%Setting model parameters:
h           = 1e-2;     %regularization function parameter  []
d60         = par{1};   %Oil specific gravity 				[]
Tbn         = par{2};   %Oil normal boiling point 			[K]
Cws     	= par{3};	%Water salt concentration			[g/L]	
dg          = par{4};   %Gas specific gravity 				[]
Cv          = par{5};   %Valve constant US units 			[gpm]
FL          = par{6};   %Liquid pressure recovery factor 	[]

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

%Valve coefficient constant [m^2]
Av=Cv*(0.0000630902)*sqrt(PP.rho_H2Ostd*14.696/101325);

switch flag,
%----- S-FUNCTION BUILDING PARAMETERS -------------------------------
    case 0,
sys(1)=0;%Number of continuous states
sys(2)=0;%Number of discrete states
sys(3)=4;%Number of outputs
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
%----- State derivatives for sys output flag=1
sys=[];
%--------------------------------------------------------------------
    case 2,
sys=[];
%----- OUPUTS -------------------------------------------------------
    case 3,        
%----- Inputs:
BSW_in	= u(1);
P_in    = u(2);
T_in	= u(3);
P_out   = u(4);
c       = u(5);

%----- Model internal variables:
%Inlet solution gas/oil ratio [Sm^3/Sm^3]
Rs_in=SolGasOilRatio(PP,P_in,T_in,d60,dg);

%Inlet oil formation volume factor [m^3/Sm^3]
Bo_in=OilFVF(PP,P_in,T_in,d60,dg);

%Inlet water formation volume factor [m^3/Sm^3]
Bw_in=WaterFVF(PP,P_in,T_in);

%Inlet liquid phase mass density [kg/m^3]
rL_in=((1-BSW_in)*d60*PP.rho_H2Ostd+BSW_in*rWstd)/...
 ((1-BSW_in)*Bo_in+BSW_in*Bw_in);

%Pressure drop [Kpa]
DP=P_in-P_out;

%Liquid phase vapour pressure [kPa]
PvL=P_in;

%Oil phase pseudo-critical pressure [kPa]
%RIAZI (2005) eq. (2.68) pp. 60
PpcO=(100)*6.9575*exp(-1.35e-2*Tbn-0.3129*d60+9.174e-3*Tbn*d60)*...
 (Tbn)^0.6791/d60^0.6807;

%Oil phase molar mass [kg/kmol]
%API-TDB Procedure 2B2.1 (Daubert and Danner 1997)
MwO=20.486*(Tbn*1.8)^1.26007*d60^4.98308*...
 exp(1.165e-4*(Tbn*1.8)-7.78712*d60+1.1582e-3*(Tbn*1.8)*d60);

%Molar mixing rule for liquid phase pseudo-critical pressure [kPa]
PpcL=(PpcO*(1-BSW_in)*d60*PP.rho_H2Ostd/MwO+22064*BSW_in*rWstd/18.02)/...
 ((1-BSW_in)*d60*PP.rho_H2Ostd/MwO+BSW_in*rWstd/18.02);

%Linear intrinsic characteristic
fc=c;

%Liquid critical pressure ratio factor
%{
FF=0.5*(1+tanh((PvL-PpcL)/(100)/1e-6))*(0.96-0.28) +...
 0.5*(1+tanh((PpcL-PvL)/(100)/1e-6))*(0.96-0.28*sqrt(PvL/PpcL));
%}
FF=0.5*(1+tanh((PvL-PpcL)/(100)/h))*(0.96-0.28) +...
 0.5*(1+tanh((PpcL-PvL)/(100)/h))*(0.96-0.28*sqrt(PvL/PpcL));

%Pressure drop at choked flow
DPchok=FL^2*(P_in-FF*PvL);

%----- Output equations:
%Pressure-flow relation
%{
DPflow=0.5*(1+tanh(((DPchok-DP)+1e-3)/1e-4))*DP...
 +0.5*(1+tanh((DP-DPchok)/1e-4))*DPchok;
%}
DPflow=0.5*(1+tanh(((DPchok-DP)+1e-3)/h))*DP...
 +0.5*(1+tanh((DP-DPchok)/h))*DPchok;
Qf=(24*3600)/((1-BSW_in)*d60*PP.rho_H2Ostd+BSW_in*rWstd)*...
 fc*Av*sqrt(rL_in*DPflow*1e3);

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
%Qg=0.5*(1+tanh(-(Rseq-Rs_in)/1e-6))*Qf*(1-BSW)*(Rs_in-Rseq);
Qg=0.5*(1+tanh(-(Rseq-Rs_in)/h))*Qf*(1-BSW)*(Rs_in-Rseq);

%Outlet gas specific residual enthalpy [kJ/kg]
eGres=GasSpecResEnth(PP,P_out,T,dg);

%Outlet temperature [K]
%{
T = PP.Tref + (...
	Qf*(1-BSW)*d60*PP.rho_H2Ostd*cPo*(T_in-PP.Tref) ...
	+ Qf*BSW*rWstd*cPw*(T_in-PP.Tref) ...
	- Qg*PP.rho_AIRstd*dg*(eGres+DeOG))/(...
		Qf*((1-BSW)*d60*PP.rho_H2Ostd*cPo + BSW*rWstd*cPw) + ...
		Qg*PP.rho_AIRstd*dg*cPg);
%}
T=T_in*((1-BSW)*d60*PP.rho_H2Ostd*cPo+BSW*rWstd*cPw)/...
    ((1-BSW)*d60*PP.rho_H2Ostd*cPo+BSW*rWstd*cPw+Qg/Qf*dg*PP.rho_AIRstd*cPg)...
 - Qg/Qf*dg*PP.rho_AIRstd/...
    ((1-BSW)*d60*PP.rho_H2Ostd*cPo+BSW*rWstd*cPw+Qg/Qf*dg*PP.rho_AIRstd*cPg)*...
        (eGres + DeOG - cPg*PP.Tref);

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