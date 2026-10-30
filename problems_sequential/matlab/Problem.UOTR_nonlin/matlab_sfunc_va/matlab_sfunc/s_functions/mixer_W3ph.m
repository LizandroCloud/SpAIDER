function [sys,x0,str,ts,simStateCompliance]=mixer_W3ph(t,x,u,flag,cinit,par)
%mixer_W3ph - Water x 3-phase streams mixer non-linear model S-Function
%
%   This model describes the dynamic behavior of a mixer with two ty-
%	pes ofinlet streams: (F) triphase oil-gas-water and; (W) water 
%	brine only for oil & gas upstream process, using the black-oil 
%	approach.
%
%   [sys,x0,str,ts,simStateCompliance]=mixer_W3ph(t,x,u,flag,cinit,par)
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
%   par{5}: nF      Number of F inlet streams                []
%   par{6}: nW      Number of W inlet streams                []
%
%   States:
%   There are no states for this model.
%
%   Inputs:
%u(4*(k-1)+1):  	F(k).Qf  Inlet F liquid std. vol. flow rate [Sm^3/d]
%u(4*(k-1)+2):  	F(k).Qg  Inlet F gas std. vol. flow rate    [Sm^3/d]
%u(4*(k-1)+3):  	F(k).BSW Inlet F basic sediments & water    []
%u(4*(k-1)+4):  	F(k).T   Inlet F temperature                [K]
%u(4*nF+2*(k-1)+1): W(k).Qw  Inlet W liquid std. vol. flow rate [Sm^3/d]
%u(4*nF+2*(k-1)+2): W(k).T   Inlet W temperature                [K]
%u(4*nF+2*nW+1): 	P        Mixer node pressure                [kPa]
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
d60         = par{1};   %Oil specific gravity 				[]
Tbn         = par{2};   %Oil normal boiling point 			[K]
Cws     	= par{3};	%Water salt concentration			[g/L]	
dg          = par{4};   %Gas specific gravity 				[]
nF          = par{5};   %Number of F inlet streams 			[]
nW          = par{6};   %Number of W inlet streams 			[]

%Water salt mass fraction []
ws=Cws/(PP.rho_H2Ostd+Cws);

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
sys(4)=4*nF+2*nW+1;%Number of inputs
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
for k=1:nF
F(k).Qf		= u(4*(k-1)+1);
F(k).Qg		= u(4*(k-1)+2);
F(k).BSW	= u(4*(k-1)+3);
F(k).T     	= u(4*(k-1)+4);
end
for k=1:nW
W(k).Qw		= u(4*nF+2*(k-1)+1);
W(k).T     	= u(4*nF+2*(k-1)+2);
end
P       	= u(4*nF+2*nW+1);

%----- Model internal variables:
for k=1:nF
%Inlet F solution gas/oil ratio [Sm^3/Sm^3]
F(k).Rs=SolGasOilRatio(PP,P,F(k).T,d60,dg);

%Inlet F gas specific residual enthalpy [kJ/kg]
F(k).eGres=GasSpecResEnth(PP,P,F(k).T,dg);
end

%----- Output equations:
%Liquid std. vol. balance
Qf=sum([F.Qf])+sum([W.Qw]);

%Basic sediments & water
BSW=(sum([F.Qf].*[F.BSW])+sum([W.Qw]))/Qf;

%Energy balance: isenthalpic expansion - Iterative calculation
T=min([[F.T],[W.T]]); iter=0; maxiter=100; rtol=1e-5;
while iter<maxiter
    iter=iter+1;
    T0=T;
    
%Solution gas/oil ratio [Sm^3/Sm^3]
Rseq=SolGasOilRatio(PP,P,T,d60,dg);

%Sum of inlet gas std. vol.
Qg_sumIn=sum([F.Qf].*(1-[F.BSW]).*[F.Rs]+[F.Qg]);

%Gas std. vol. unbalance
Qg_unbal=Qf*(1-BSW)*Rseq-Qg_sumIn;

%Std. vol. balance for liberated gas [Sm^3/d]
Qg=0.5*(1+tanh(-(Qg_unbal)/1e-6))*(-Qg_unbal);

%Outlet gas specific residual enthalpy [kJ/kg]
eGres=GasSpecResEnth(PP,P,T,dg);

%Outlet temperature [K]
%{
%Gas specific residual enthalpy derivative do T [kJ/kg K]
deGresdT=dGasSpecResEnthdT(PP,P,T,dg);

f=Qf*((1-BSW)*d60*PP.rho_H2Ostd*cPo+BSW*rWstd*cPw)*(T-PP.Tref)...
 + Qg*PP.rho_AIRstd*dg*(eGres+DeOG+cPg*(T-PP.Tref)) ...
 - sum([F.Qf].*((1-[F.BSW])*d60*PP.rho_H2Ostd*cPo+[F.BSW]*rWstd*cPw).*([F.T]-PP.Tref))...
 - sum([F.Qg]*PP.rho_AIRstd*dg.*([F.eGres]+DeOG+cPg.*([F.T]-PP.Tref)))...
 - sum([W.Qw]*rWstd*cPw.*([W.T]-PP.Tref));
dfdT=Qf*((1-BSW)*d60*PP.rho_H2Ostd*cPo+BSW*rWstd*cPw)...
 + Qg*PP.rho_AIRstd*dg*(deGresdT+cPg);
T=T-f/dfdT;
%}

T = PP.Tref + (...
	sum([F.Qf].*(1-[F.BSW])*d60*PP.rho_H2Ostd*cPo.*([F.T]-PP.Tref)) ...
	+ sum([F.Qf].*[F.BSW]*rWstd*cPw.*([F.T]-PP.Tref)) ...
	+ sum([F.Qg]*PP.rho_AIRstd*dg.*([F.eGres]+DeOG+cPg.*([F.T]-PP.Tref))) ...
	+ sum([W.Qw]*rWstd*cPw.*([W.T]-PP.Tref)) ...
	- Qg*PP.rho_AIRstd*dg*(eGres+DeOG))/(...
		Qf*((1-BSW)*d60*PP.rho_H2Ostd*cPo + BSW*rWstd*cPw) ...
		+ Qg*PP.rho_AIRstd*dg*cPg);

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