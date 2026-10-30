function [sys,x0,str,ts,simStateCompliance]=oil_treater(t,x,u,flag,cinit,par)
%oil_treater - oil treater vessel non-linear model S-Function
%
%   This model describes the dynamic behavior of an upstream oil &
%   gas process oil treater vessel for oil-water separation, using 
%   the black-oil approach.
%
%   [sys,x0,str,ts,simStateCompliance]=oil_treater(t,x,u,flag,cinit,par)
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
%   par{5}:  C       Vessel length                            	[m]
%   par{6}:  D       Vessel diameter                          	[m]
%   par{7}:  tRL50   Res. time for water half-removal    		[s]
%   par{8}:  gamma   Water removal eff. correl. 2nd par 		[]
%
%   Canonical states:
%   x(1): VWstd     Water phase std. vol. hold-up            	[m^3]
%   x(2): E         Energy hold-up                           	[kJ]
%
%   Process variable states:
%   xpv(1): hW      Water phase level                        	[m]
%   xpv(2): T       Temperature                              	[K]
%
%   Inputs:
%   u(1): BSW_in    Inlet basic sediments & water            []
%   u(2): P_in      Inlet pressure                           [kPa]
%   u(3): T_in      Inlet temperature                        [K]
%   u(4): Qf        Light phase liq std. vol. out flow rate  [Sm^3/d]
%   u(5): Qw        Water phase liq std. vol. out flow rate  [Sm^3/d]
%
%   Outputs:
%   y(1): hW        Water phase level                        [m]
%   y(2): BSW       Basic sediments & water                  []
%   y(3): Qf_in     Inlet std. vol. flow rate                [Sm^3/d]
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
tRL50   = par{7};   %Residence time for water half-removal 	[s]
gamma   = par{8};   %water removal eff. correl. 2nd parameter []

%Water salt mass fraction []
ws=Cws/(PP.rho_H2Ostd+Cws);

%Water brine mass density [kg/m^3]
rWstd=WaterStdMassDens(PP,ws);

%Oil specific heat capacity [kJ/kg K]
cPo=OilSpecHeatCap(PP,d60,Tbn);

%Water specific heat capacity [kJ/kg K]
cPw=WaterSpecHeatCap(PP,PP.Pref,PP.Tref,ws);
   
%Vessel total volume [m^3]   
Vt=0.25*C*pi*D^2;
   
switch flag,
%----- S-FUNCTION BUILDING PARAMETERS AND INITIAL CONDITIONS --------
    case 0,        
sys(1)=2;%Number of continuous states
sys(2)=0;%Number of discrete states
sys(3)=4;%Number of outputs
sys(4)=5;%Number of inputs
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
hW 		= x(1);
T		= x(2);

%----- Inputs:
BSW_in	= u(1);
P_in	= u(2);
T_in	= u(3);
Qf      = u(4);
Qw      = u(5);

%----- Model internal variables:
%Oil formation volume factor [m^3/Sm^3]
Bo=OilFVF(PP,P_in,T,d60,dg);

%Water formation volume factor [m^3/Sm^3]
Bw=WaterFVF(PP,P_in,T);

%Water brine FVF derivative to T
dBwdT=dWaterFVFdT(PP,P_in,T);

%Swetted cross section area [m^2]
Aw=D^2/8*(2*acos((D-2*hW)/D)-sin(2*acos((D-2*hW)/D)));

%BSW [] and Qf_in [Sm^3/d]
y=oil_treater(t,x,u,3,[],par);
BSW 	= y(2);
Qf_in 	= y(3);

%Water removal efficiency from light phase
effWL=(Qf_in*BSW_in-Qf*BSW)/(Qf_in*BSW_in);

%----- Canonical states:
%Water phase std. vol. hold-up [Sm^3]
VWstd=C*Aw/Bw;

%Light phase std. vol. hold-ups [Sm^3]
VLstd=(Vt-C*Aw)/((1-BSW)*Bo+BSW*Bw);
VOLstd=VLstd*(1-BSW);
VWLstd=VLstd*BSW;

%----- Canonical state equations:
%Oil std. vol. balance in liquid phase
dVOLstddt=(Qf_in*(1-BSW_in)-Qf*(1-BSW))/(24*3600);

%Water std. vol. balance in liquid phase
dVWLstddt=(Qf_in*BSW_in*(1-effWL)-Qf*BSW)/(24*3600);

%Water phase std. vol. balance
dVWstddt=(Qf_in*BSW_in*effWL-Qw)/(24*3600);

%Energy balance
dEdt=...
	( Qf_in*(1-BSW_in)*d60*PP.rho_H2Ostd*cPo*(T_in-PP.Tref)...
	+ Qf_in*BSW_in*rWstd*cPw*(T_in-PP.Tref)...
	- Qf*(1-BSW)*d60*PP.rho_H2Ostd*cPo*(T-PP.Tref)...
	- Qf*BSW*rWstd*cPw*(T-PP.Tref)...
	- Qw*rWstd*cPw*(T-PP.Tref))/(24*3600);

%----- Process variable state equations:
%2x2 linear system:
%vec=[dhWdt; dTdt];
M=[2*C*sqrt(hW*(D-hW)) 	-VWstd*dBwdT
   0  					(VWstd+VWLstd)*rWstd*cPw+VOLstd*d60*PP.rho_H2Ostd*cPo];
b=[Bw*dVWstddt; ...
   dEdt-rWstd*cPw*(T-PP.Tref)*(dVWstddt+dVWLstddt)-d60*PP.rho_H2Ostd*cPo*(T-PP.Tref)*dVOLstddt];
%vec=M\b;
vec=linsolve(M,b);
dhWdt=vec(1);
dTdt=vec(2);

%----- State derivatives for sys output flag=1
sys(1)=dhWdt;
sys(2)=dTdt;
sys=sys(:);
%--------------------------------------------------------------------
    case 2,
sys=[];
%----- OUPUTS -------------------------------------------------------
    case 3,
%----- States:
hW 		= x(1);
T		= x(2);

%----- Inputs:
BSW_in	= u(1);
P_in	= u(2);
T_in	= u(3);
Qf      = u(4);
Qw      = u(5);

%----- Model internal variables:
%Inlet oil formation volume factor [m^3/Sm^3]
Bo_in=OilFVF(PP,P_in,T_in,d60,dg);

%Oil formation volume factor [m^3/Sm^3]
Bo=OilFVF(PP,P_in,T,d60,dg);

%Water formation volume factor [m^3/Sm^3]
Bw=WaterFVF(PP,P_in,T);

%Oil FVF derivative to T
dBodT=dOilFVFdT(PP,P_in,T,d60,dg);

%Water brine FVF derivative to T
dBwdT=dWaterFVFdT(PP,P_in,T);

%Swetted cross section area [m^2]
Aw=D^2/8*(2*acos((D-2*hW)/D)-sin(2*acos((D-2*hW)/D)));

%Water phase std. vol. hold-up [Sm^3]
VWstd=C*Aw/Bw;

%Iterative calculation for y=[BSW; Qf_in]
BSW=BSW_in/(1+10^gamma);
Qf_in=Qf;
iter=0; maxiter=100; rtol=1e-5;
while iter<maxiter
    iter=iter+1;
    y0=[BSW; Qf_in];

	%Light phase std. vol. hold-ups [Sm^3]
	VLstd=(Vt-C*Aw)/((1-BSW)*Bo+BSW*Bw);
	VOLstd=VLstd*(1-BSW);
	VWLstd=VLstd*BSW;

	%Oil residence time in light phase
	tRL=(VOLstd*Bo)/(Qf_in/(24*3600)*(1-BSW_in)*Bo_in);

	%Relation between BSW and residence time for light phase
	BSW=BSW_in/(1+(tRL/tRL50)^gamma);

	%Energy hold-up differentiation to T
	dEdT=VOLstd*d60*PP.rho_H2Ostd*cPo+(VWLstd+VWstd)*rWstd*cPw;

	%Energy balance
	dEdt=...
		( Qf_in*(1-BSW_in)*d60*PP.rho_H2Ostd*cPo*(T_in-PP.Tref)...
		+ Qf_in*BSW_in*rWstd*cPw*(T_in-PP.Tref)...
		- Qf*(1-BSW)*d60*PP.rho_H2Ostd*cPo*(T-PP.Tref)...
		- Qf*BSW*rWstd*cPw*(T-PP.Tref)...
		- Qw*rWstd*cPw*(T-PP.Tref))/(24*3600);

	%Inlet liquid std. vol. flow rate [Sm^3/d]
	Qf_in=(Qf*((1-BSW)*Bo+BSW*Bw)+Qw*Bw...
		-(VOLstd*dBodT+(VWLstd+VWstd)*dBwdT)/dEdT*dEdt*(24*3600))...
		/((1-BSW_in)*Bo+BSW_in*Bw);
    
	y=[BSW; Qf_in]; rdev=abs((y-y0)./y0);
    if rdev<rtol
        break
    end
end

%----- Outputs for sys output flag=3:
sys(1)=hW;
sys(2)=BSW;
sys(3)=Qf_in;
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