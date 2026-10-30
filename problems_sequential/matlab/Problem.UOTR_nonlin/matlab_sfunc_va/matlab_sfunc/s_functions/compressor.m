function [sys,x0,str,ts,simStateCompliance]=compressor(t,x,u,flag,cinit,par)
%compressor - Gas compressor non-linear model S-Function
%
%   This model describes the dynamic behavior of a compressor with a 
%	gas-only inlet/outlet streams for oil & gas upstream processes,
%	using the black-oil approach.
%
%   [sys,x0,str,ts,simStateCompliance]=compressor(t,x,u,flag,cinit,par)
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
%   par{1}: dg      	Gas specific gravity                []
%   par{2}: Fgm     	Specific heat capacity ratio factor []
%   par{3}: Nr_Max    	Maximum compressor rotation         [rpm]
%   par{4}: PrCC  		Pres. ratio charac. curve par.   	[]
%   par{5}: EffSPOL		Isentropic eff. polynomial par.   	[]
%
%   States:
%   There are no states for this model.
%
%   Inputs:
%   u(1): P_in      	Inlet pressure                      [kPa]
%   u(2): T_in      	Inlet temperature                   [K]
%   u(3): P_out     	Outlet pressure                     [kPa]
%   u(4): MV_in     	Input signal from controller        []
%
%   Outputs:
%   y(1): Qg        	Gas std. vol. flow rate             [Sm^3/d]
%   y(2): T         	Temperature                         [K]
%	y(3): W				Compressor power 				 	[kW]
%
%   Initial conditions set:
%   There are no initial conditions for this model.

%   2023. Version 0. Coded by Diego Trica (diegotrica@gmail.com)

%Blackoil object class property package. On the main script the
%user must declare PP as global and PP=blackoil.
global PP

%Setting model parameters:
dg          = par{1};   %Gas specific gravity []
Fgm         = par{2};   %Specific heat capacity ratio factor []
Nr_Max      = par{3};   %Maximum compressor rotation [rpm]
PrCC       	= par{4};   %Pres. ratio charac. curve par. []
EffSPOL     = par{5};   %Isentropic eff. polynomial par. []

%Gas specific heat capacity [kJ/kg K]
cPg=GasSpecHeatCap(PP,dg);

switch flag,
%----- S-FUNCTION BUILDING PARAMETERS -------------------------------
    case 0,
sys(1)=0;   %Number of continuous states
sys(2)=0;   %Number of discrete states
sys(3)=3;   %Number of outputs
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
MV_in   = u(4);

%----- Model internal variables:
%Compressor speed [rpm]
Nr = MV_in*Nr_Max;

%Inlet gas compressibility factor []
Zg_in=GasCompFactor(PP,P_in,T_in,dg);

%Inlet gas mass density [kg/m^3]
%rG_in=P_in*(dg*PP.MM_Air)/(PP.R*Zg_in*T_in);

%Inlet gas specific residual enthalpy [kJ/kg]
eGres_in=GasSpecResEnth(PP,P_in,T_in,dg);

%Pressure ratio []
Pr=P_out/P_in;

%----- Output equations:
%Pressure-flow relation - CAMPOS & TEIXEIRA (2010)
%Pr^((Fgm-1)/Fgm)-1=(PrCC(1)*1e-8)*(Nr)^2...
% +(PrCC(2)*1e-6)*(Nr)*(Qg/(24*3600)*dg*PP.rho_AIRstd)/(P_in*1e-2)...
% -(PrCC(3)*1e-1)*(Qg/(24*3600)*dg*PP.rho_AIRstd/(P_in*1e-2))^2;
a=-PrCC(3)*1e-1;
b=PrCC(2)*1e-6*Nr;
c=1-Pr^((Fgm-1)/Fgm)+PrCC(1)*1e-8*Nr^2;
sol=-(b/(2*a))+sqrt((b/(2*a))^2-c/a); %polynomial solution
Mg=sol*(P_in*1e-2);  				%Mass flow rate [kg/s]
Qg=(24*3600)*Mg/(dg*PP.rho_AIRstd); %Std. vol. flow rate [Sm^3/d]

%Mass flow rate x isentropic efficiency: Quadratic polynomial
effS=EffSPOL(1)+(EffSPOL(2)*1e-1)*Mg-(EffSPOL(3)*1e-3)*Mg^2;

%Energy balance: isentropic compression - Iterative calculation
T=T_in+50; iter=0; maxiter=100; rtol=1e-5;
while iter<maxiter
    iter=iter+1;
    T0=T;
	
%Gas compressibility factor []
Zg=GasCompFactor(PP,P_out,T,dg);	

%Isentropic head [kJ/kg K]
HeadS=Fgm/(Fgm-1)*PP.R/(dg*PP.MM_Air)*...
 T_in*0.5*(Zg_in+Zg)*(Pr^((Fgm-1)/Fgm)-1);
 
%Outlet gas specific residual enthalpy [kJ/kg]
eGres=GasSpecResEnth(PP,P_out,T,dg);

%Energy balance: Outlet temperature [K]
T=T_in+(eGres_in+(HeadS/effS)-eGres)/cPg;

	rdev=abs((T-T0)/T0);
	if rdev<rtol
		break
	end
end

%Compressor fluid power [kW]
W=(HeadS/effS)*Qg/(24*3600)*dg*PP.rho_AIRstd;

%----- Outputs for sys output flag=3:
sys(1)=Qg;
sys(2)=T;
sys(3)=W;
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