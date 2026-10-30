function [sys,x0,str,ts,simStateCompliance]=degasser(t,x,u,flag,cinit,par)
%degasser - degasser vessel non-linear model S-Function
%
%   This model describes the dynamic behavior of an upstream oil &
%   gas process degasser vessel for gas-liquid separation, using the 
%   black-oil approach.
%
%   [sys,x0,str,ts,simStateCompliance]=degasser(t,x,u,flag,cinit,par)
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
%   par{1}:  d60     Oil specific gravity                     	[]
%   par{2}:  Tbn     Oil normal boiling point                 	[K]
%   par{3}:  Cws     Water salt concentration                 	[g/L]
%   par{4}:  dg      Gas specific gravity                     	[]
%   par{5}:  C       Vessel length                             	[m]
%   par{6}:  D       Vessel diameter                          	[m]
%
%   Canonical states:
%   x(1): VOLstd    Oil std. vol. hold-up in liquid phase    	[m^3]
%   x(2): VWLstd    Water std. vol. hold-up in liquid phase  	[m^3]
%   x(3): VGstd     Gas std. vol. hold-up                    	[m^3]
%   x(4): E         Energy hold-up                           	[kJ]
%
%   Process variable states:
%   xpv(1): hL      Liquid level                       			[m]
%   xpv(2): BSW	 	Basic sediments & water 					[]
%   xpv(3): P       Pressure                                 	[kPa]
%   xpv(4): T       Temperature                              	[K]
%
%   Inputs:
%   u(1): Qf_in     Inlet liquid std. vol. flow rate         [Sm^3/d]
%   u(2): Qg_in     Inlet gas std. vol. flow rate            [Sm^3/d]
%   u(3): BSW_in    Inlet basic sediments & water            []
%   u(4): T_in      Inlet temperature                        [K]
%   u(5): Qf        Liquid std. vol. flow rate               [Sm^3/d]
%   u(6): Qg        Gas std. vol. flow rate                  [Sm^3/d]
%
%   Outputs:
%   y(1): hL        Liquid level                             [m]
%   y(2): BSW       Basic sediments & water                  []
%   y(3): P         Pressure                                 [kPa]
%   y(4): T         Temperature                              [K]
%
%   Initial conditions set:
%   cinit = xpv

%   2022. Coded by Diego Trica (diegotrica@gmail.com)

%Blackoil object class property package. On the main script the
%user must declare PP as global and PP=blackoil.
global PP

%Setting model parameters:
d60     = par{1};	%Oil specific gravity 					[]
Tbn     = par{2};	%Oil normal boiling point 				[K]
Cws     = par{3};	%Water salt concentration				[g/L]
dg      = par{4};	%Gas specific gravity 					[]    
C       = par{5};   %Vessel length 							[m]
D       = par{6};   %Vessel diameter 						[m]

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
   
%Vessel total volume [m^3]   
Vt=0.25*C*pi*D^2;
   
switch flag,
%----- S-FUNCTION BUILDING PARAMETERS AND INITIAL CONDITIONS --------
    case 0,        
sys(1)=4;%Number of continuous states
sys(2)=0;%Number of discrete states
sys(3)=4;%Number of outputs
sys(4)=6;%Number of inputs
sys(5)=0;%Reserved for root finding. Must be zero.
sys(6)=1;%Direct Feedthrough S-function parameter
sys(7)=1;%Number of sample times (at least one sample time is needed)

%----- States initial condition:
x0=cinit;
x0=x0(:);

%----- Other building parameters:
str=[];     %str is always an empty matrix
ts=[0 0];   %Sample time array [0 0] for continuous system
simStateCompliance='UnknownSimState';
%--------------------------------------------------------------------
%----- STATE DERIVATIVES --------------------------------------------
    case 1,
%----- States:
hL  	= x(1);
BSW 	= x(2);
P       = x(3);
T		= x(4);

%----- Inputs:
Qf_in	= u(1);
Qg_in	= u(2);
BSW_in	= u(3);
T_in	= u(4);
Qf      = u(5);
Qg      = u(6);

%----- Model internal variables:
%Oil formation volume factor [m^3/Sm^3]
Bo=OilFVF(PP,P,T,d60,dg);

%Water formation volume factor [m^3/Sm^3]
Bw=WaterFVF(PP,P,T);

%Gas formation volume factor [m^3/Sm^3]
Bg=GasFVF(PP,P,T,dg);

%Inlet solution gas/oil ratio [Sm^3/Sm^3]
Rs_in=SolGasOilRatio(PP,P,T_in,d60,dg);

%Solution gas-oil ratio [Sm^3/Sm^3]
Rs=SolGasOilRatio(PP,P,T,d60,dg);

%Inlet gas specific residual enthalpy [kJ/kg]
eGres_in=GasSpecResEnth(PP,P,T_in,dg);

%Outlet gas specific residual enthalpy [kJ/kg]
eGres=GasSpecResEnth(PP,P,T,dg);

%Oil FVF derivative to P and T
dBodP=dOilFVFdP(PP,P,T,d60,dg);
dBodT=dOilFVFdT(PP,P,T,d60,dg);

%Water brine FVF derivative to P and T
dBwdP=dWaterFVFdP(PP,P,T);
dBwdT=dWaterFVFdT(PP,P,T);

%Gas FVF derivative to P and T
dBgdP=dGasFVFdP(PP,P,T,dg);
dBgdT=dGasFVFdT(PP,P,T,dg);

%Swetted cross section area [m^2]
Aliq=D^2/8*(2*acos((D-2*hL)/D)-sin(2*acos((D-2*hL)/D)));

%----- Canonical states:
%2x2 linear system:
%vec=[VOLstd; VWLstd];
%(1) VOLstd*Bo+VWLstd*Bw=C*Aliq;
%(2) BSW*VOLstd+(BSW-1)*VWLstd=0;
M=[Bo  Bw
   BSW (BSW-1)];
b=[C*Aliq; 0];
vec=M\b;
VOLstd=vec(1);
VWLstd=vec(2);

%Gas phase std. vol. hold-up [Sm^3]
VGstd=(Vt-C*Aliq)/Bg;

%----- Canonical state equations:
%Oil std. vol. balance
dVOLstddt=(Qf_in*(1-BSW_in)-Qf*(1-BSW))/(24*3600);

%Water std. vol. balance
dVWLstddt=(Qf_in*BSW_in-Qf*BSW)/(24*3600);

%Gas std. vol. balance
dVGstddt=(Qf_in*(1-BSW_in)*Rs_in+Qg_in-Qf*(1-BSW)*Rs-Qg)/(24*3600);

%Energy balance
dEdt=...
 ( Qf_in*(1-BSW_in)*d60*PP.rho_H2Ostd*cPo*(T_in-PP.Tref)...
 + Qf_in*BSW_in*rWstd*cPw*(T_in-PP.Tref)...
 + Qg_in*dg*PP.rho_AIRstd*(eGres_in + DeOG + cPg*(T_in-PP.Tref))...
 - Qf*(1-BSW)*d60*PP.rho_H2Ostd*cPo*(T-PP.Tref)...
 - Qf*BSW*rWstd*cPw*(T-PP.Tref)...
 - Qg*dg*PP.rho_AIRstd*(eGres + DeOG + cPg*(T-PP.Tref)))/(24*3600);

%----- Process variable state equations:
%3x3 linear system:
%vec=[dhLdt; dPdt; dTdt];
M=[2*C*sqrt(hL*(D-hL)) 	-VWLstd*dBwdP-VOLstd*dBodP  				-VWLstd*dBwdT-VOLstd*dBodT
   0  					VGstd*dBgdP+VWLstd*dBwdP+VOLstd*dBodP   	VGstd*dBgdT+VWLstd*dBwdT+VOLstd*dBodT
   0  					0 											VWLstd*rWstd*cPw+VOLstd*d60*PP.rho_H2Ostd*cPo];
b=[Bw*dVWLstddt+Bo*dVOLstddt; ...
   -Bg*dVGstddt-Bw*dVWLstddt-Bo*dVOLstddt; ...
   dEdt-rWstd*cPw*(T-PP.Tref)*dVWLstddt-d60*PP.rho_H2Ostd*cPo*(T-PP.Tref)*dVOLstddt];
vec=M\b;
dhLdt=vec(1);
dPdt=vec(2);
dTdt=vec(3);

%Water std. vol. content at light phase in CS rate:
dBSWdt=((1-BSW)*dVWLstddt-BSW*dVOLstddt)/(VOLstd+VWLstd);

%----- State derivatives for sys output flag=1
sys(1)=dhLdt;
sys(2)=dBSWdt;
sys(3)=dPdt;
sys(4)=dTdt;
sys=sys(:);
%--------------------------------------------------------------------
    case 2,
sys=[];
%----- OUPUTS -------------------------------------------------------
    case 3,
%----- Outputs for sys output flag=3:
sys=x;
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