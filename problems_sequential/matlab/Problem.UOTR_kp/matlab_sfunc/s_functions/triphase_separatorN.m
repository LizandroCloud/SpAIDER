function [sys,x0,str,ts,simStateCompliance]=triphase_separator(t,x,u,flag,cinit,par)
%triphase_separator - triphase separator vessel non-linear model
%S-Function
%
%   This model describes the dynamic behavior of an upstream oil &
%   gas process triphase separator vessel, using the black-oil approach.
%
%   [sys,x0,str,ts,simStateCompliance]=triphase_separator(t,x,u,flag,cinit,par)
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
%   par{5}:  LCS     Separation chamber (CS) length           	[m]
%   par{6}:  LCL     Lights chamber (CL) length               	[m]
%   par{7}:  D       Vessel diameter                          	[m]
%   par{8}:  hweir   Weir height               				 	[m]
%
%   Canonical states:
%   xcs(1): VOLCSstd  Oil std. vol. hold-up in CS liq. phase   [Sm^3]
%   xcs(2): VWLCSstd  Water std. vol. hold-up in CS liq. phase [Sm^3]
%   xcs(3): VWCSstd   Water phase std. vol. hold-up in CS      [Sm^3]
%   xcs(4): VOLCLstd  Oil std. vol. hold-up in CL liq. phase   [Sm^3]
%   xcs(5): VWLCLstd  Water std. vol. hold-up in CL liq. phase [Sm^3]
%   xcs(6): VGstd     Gas std. vol. hold-up                    [Sm^3]
%   xcs(7): E         Energy hold-up                           [kJ]
%
%   Process variable states:
%   xpv(1): hWCS      Water phase level in CS                  	[m]
%   xpv(2): BSWLCS    Basic sediments & water in CS 			[]
%   xpv(3): hCS       Total level in CS                  		[m]
%   xpv(4): hCL       Liquid level in CL                       	[m]
%   xpv(5): BSWLCL	  Basic sediments & water in CL 			[]
%   xpv(6): P         Pressure                                 	[kPa]
%   xpv(7): T         Temperature                              	[K]
%
%   Inputs:
%   u(1): Qf_in     Inlet liquid std. vol. flow rate         [Sm^3/d]
%   u(2): Qg_in     Inlet gas std. vol. flow rate            [Sm^3/d]
%   u(3): BSW_in    Inlet basic sediments & water            []
%   u(4): T_in      Inlet temperature                        [K]
%   u(5): Qf        Liquid std. vol. flow rate               [Sm^3/d]
%   u(6): Qw        Water std. vol. flow rate                [Sm^3/d]
%   u(7): Qg        Gas std. vol. flow rate                  [Sm^3/d]
%
%   Outputs:
%   y(1): hWCS      Water phase level in CS   					[m]
%   y(2): BSWLCS    Basic sediments & water in CS 				[]
%   y(3): hCS       Total liquid level in CS                  	[m]
%   y(4): hCL       Liquid level in CL                       	[m]
%   y(5): BSWLCL	Basic sediments & water in CL 				[]
%   y(6): P         Pressure                                 	[kPa]
%   y(7): T         Temperature                              	[K]
%
%   Initial conditions set:
%   cinit = xpv

%   2023. Coded by Diego Trica (diegotrica@gmail.com)

%Blackoil object class property package. On the main script the
%user must declare PP as global and PP=blackoil.
global PP

%----- Parameters:
h       = 1e-2;     %regularization function parameter  []
d60     = par{1};	%Oil specific gravity 					[]
Tbn     = par{2};	%Oil normal boiling point 				[K]
Cws     = par{3};	%Water salt concentration				[g/L]
dg      = par{4};	%Gas specific gravity 					[]
LCS     = par{5};	%Separation chamber (CS) length 		[m]
LCL     = par{6};	%Lights chamber (CL) length 			[m]
D       = par{7};	%Vessel diameter 						[m]
hweir   = par{8};	%Weir height 							[m]
	
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
Vt=0.25*(LCS+LCL)*pi*D^2;

%Weir length [m]
Lweir=2*sqrt(hweir*(D-hweir));

%----- States:
hWCS  	= x(1);
BSWLCS 	= x(2);
hCS  	= x(3);
hCL  	= x(4);
BSWLCL 	= x(5);
P       = x(6);
T		= x(7);

%----- Inputs:
Qf_in	= u(1);
Qg_in	= u(2);
BSW_in	= u(3);
T_in	= u(4);
Qf      = u(5);
Qw      = u(6);
Qg      = u(7);

if flag==1 || flag==10 || flag==11
%----- Internal calculations:
%Inlet oil formation volume factor [m^3/Sm^3]
Bo_in=OilFVF(PP,P,T_in,d60,dg);

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

%Weir std. vol. flow rate [Sm^3/d]
%{
Qweir=(24*3600)/((1-BSWLCS)*Bo+BSWLCS*Bw)*( ...
 0.5*(1+tanh((hweir-hCL)/1e-4))*0.5*(1+tanh((hCS-hweir)/1e-4))* ...
  (24.88*sqrt(2*9.81)*(Lweir-0.2*(abs(hCS-hweir)))*(abs(hCS-hweir))^1.5) ...
 +0.5*(1+tanh((hCL-hweir)/1e-4))*0.5*(1+tanh((hCS-hCL)/1e-4))* ...
  (24.88*sqrt(2*9.81)*(2*sqrt(hCL*(D-hCL))-0.2*(abs(hCS-hCL)))*(abs(hCS-hCL))^1.5) ...
 );
%}
Qweir=(24*3600)/((1-BSWLCS)*Bo+BSWLCS*Bw)*( ...
 0.5*(1+tanh((hweir-hCL)/h))*0.5*(1+tanh((hCS-hweir)/h))* ...
  (24.88*sqrt(2*9.81)*(Lweir-0.2*(abs(hCS-hweir)))*(abs(hCS-hweir))^1.5) ...
 +0.5*(1+tanh((hCL-hweir)/h))*0.5*(1+tanh((hCS-hCL)/h))* ...
  (24.88*sqrt(2*9.81)*(2*sqrt(hCL*(D-hCL))-0.2*(abs(hCS-hCL)))*(abs(hCS-hCL))^1.5) ...
 );

%Water removal efficiency from light phase
effWL=0.5486-0.2675*(hWCS)-0.3626*(Qf_in/(24*3600)*(1-BSW_in)*Bo_in)...
 +0.1356*(hWCS)*(Qf_in/(24*3600)*(1-BSW_in)*Bo_in)...
 -0.0747*(hWCS)^2+0.089*(Qf_in/(24*3600)*(1-BSW_in)*Bo_in)^2;

%Water phase swetted cross section area in CS [m^2]
AWCS=D^2/8*(2*acos((D-2*hWCS)/D)-sin(2*acos((D-2*hWCS)/D)));

%Total swetted cross section area in CS [m^2]
AtCS=D^2/8*(2*acos((D-2*hCS)/D)-sin(2*acos((D-2*hCS)/D)));

%Liquid phase swetted cross section area in CL [m^2]
ALCL=D^2/8*(2*acos((D-2*hCL)/D)-sin(2*acos((D-2*hCL)/D)));

%----- Canonical states:
%5x5 linear system:
%vec=[VOLCSstd; VWLCSstd; VWCSstd; VOLCLstd; VWLCLstd];
%(1) Bw*VWCSstd=LCS*AWCS;
%(2) Bo*VOLCSstd+Bw*VWLCSstd+Bw*VWCSstd=LCS*AtCS;
%(3) Bo*VOLCLstd+Bw*VWLCLstd=LCL*ALCL;
%(4) BSWLCS*VOLCSstd+(BSWLCS-1)*VWLCSstd=0;
%(5) BSWLCL*VOLCLstd+(BSWLCL-1)*VWLCLstd=0;
H=[0   		0       	Bw 0   		0
   Bo  		Bw 	   		Bw 0   		0
   0   		0       	0  Bo  		Bw
   BSWLCS 	(BSWLCS-1) 	0  0   		0
   0   		0       	0  BSWLCL 	(BSWLCL-1)];
b=[LCS*AWCS; LCS*AtCS; LCL*ALCL; 0; 0];
%vec=H\b;
vec=linsolve(H,b);
VOLCSstd=vec(1);
VWLCSstd=vec(2);
VWCSstd=vec(3);
VOLCLstd=vec(4);
VWLCLstd=vec(5);

%Gas phase std. vol. hold-up [Sm^3]
VGstd=(Vt-(VWCSstd+VWLCSstd+VWLCLstd)*Bw-(VOLCSstd+VOLCLstd)*Bo)/Bg;

%----- Canonical state derivatives:
%Oil std. vol. balance in CS liquid phase
dVOLCSstddt=(Qf_in*(1-BSW_in)-Qweir*(1-BSWLCS))/(24*3600);

%Water std. vol. balance in CS liquid phase
dVWLCSstddt=(Qf_in*BSW_in*(1-effWL)-Qweir*BSWLCS)/(24*3600);

%Water phase std. vol. balance in CS
dVWCSstddt=(Qf_in*BSW_in*effWL-Qw)/(24*3600);

%Oil std. vol. balance in CL liquid phase
dVOLCLstddt=(Qweir*(1-BSWLCS)-Qf*(1-BSWLCL))/(24*3600);

%Water std. vol. balance in CL liquid phase
dVWLCLstddt=(Qweir*BSWLCS-Qf*BSWLCL)/(24*3600);

%Gas std. vol. balance
dVGstddt=(Qf_in*(1-BSW_in)*Rs_in+Qg_in-Qf*(1-BSWLCL)*Rs-Qg)/(24*3600);

%Energy balance
dEdt=...
 (Qf_in*(1-BSW_in)*d60*PP.rho_H2Ostd*cPo*(T_in-PP.Tref)...
 +Qf_in*BSW_in*rWstd*cPw*(T_in-PP.Tref)...
 +Qg_in*dg*PP.rho_AIRstd*(eGres_in+DeOG+cPg*(T_in-PP.Tref))...
 -Qf*(1-BSWLCL)*d60*PP.rho_H2Ostd*cPo*(T-PP.Tref)...
 -Qf*BSWLCL*rWstd*cPw*(T-PP.Tref)...
 -Qg*dg*PP.rho_AIRstd*(eGres+DeOG+cPg*(T-PP.Tref))...
 -Qw*rWstd*cPw*(T-PP.Tref))/(24*3600);

%----- Process variable derivatives:
%h(xcs,xpv,par)=0;	%algebraic relation between xpv and xcs
%dhdt = 0 -> dhdxpv*dxpvdt + dhdxcs*dxcsdt = 0;
%dxcsdt=fcs(xcs,xpv,u,par); dxpvdt=fpv(xcs,xpv,u,par);
%dhdxpv*fpv + dhdxcs*fcs = 0;
%M=dhdxpv(xcs,xpv,par); N=-dhdxcs(xcs,xpv,par);
%M*fpv - N*fcs = 0 -> fpv = M\N*fcs;
 
%dhdxpv: Algebraic relations between xcs and xpv partial derivative in relation to xpv
M=[2*LCS*sqrt(hWCS*(D-hWCS))   				0						0       				0 									0						-VWCSstd*dBwdP   															-VWCSstd*dBwdT
   0										(VOLCSstd+VWLCSstd)		0						0									0						0																			0
   -2*LCS*sqrt(hWCS*(D-hWCS))  				0						2*LCS*sqrt(hCS*(D-hCS)) 0 									0						-VWLCSstd*dBwdP-VOLCSstd*dBodP  											-VWLCSstd*dBwdT-VOLCSstd*dBodT
   0   										0						0       				2*LCL*sqrt(hCL*(D-hCL)) 			0						-VWLCLstd*dBwdP-VOLCLstd*dBodP  											-VWLCLstd*dBwdT-VOLCLstd*dBodT
   0										0						0						0									(VOLCLstd+VWLCLstd)		0																			0
   0 										0						0 						0  									0						VGstd*dBgdP+(VWCSstd+VWLCSstd+VWLCLstd)*dBwdP+(VOLCSstd+VOLCLstd)*dBodP   	VGstd*dBgdT+(VWCSstd+VWLCSstd+VWLCLstd)*dBwdT+(VOLCSstd+VOLCLstd)*dBodT
   0   										0						0       				0  									0						0 																			(VWCSstd+VWLCSstd+VWLCLstd)*rWstd*cPw+(VOLCSstd+VOLCLstd)*d60*PP.rho_H2Ostd*cPo];

%dhdxcs: Algebraic relations between xcs and xpv partial derivative in relation to xcs
N=[0 										0 						Bw 						0 									0 						0 																			0
   -BSWLCS 									(1-BSWLCS)				0						0									0						0																			0
   Bo										Bw						0						0									0						0																			0
   0										0						0						Bo									Bw						0																			0
   0 										0						0						-BSWLCL								(1-BSWLCL)				0																			0
   -Bo										-Bw						-Bw						-Bo									-Bw						-Bg																			0
   -d60*PP.rho_H2Ostd*cPo*(T-PP.Tref)		-rWstd*cPw*(T-PP.Tref)	-rWstd*cPw*(T-PP.Tref)	-d60*PP.rho_H2Ostd*cPo*(T-PP.Tref)	-rWstd*cPw*(T-PP.Tref)	0																			eye(1)]; 

%diff(xcs)
fcs=[...
	dVOLCSstddt; ...
	dVWLCSstddt; ...
	dVWCSstddt; ...
	dVOLCLstddt; ...
	dVWLCLstddt; ...
	dVGstddt; ...
	dEdt];

%diff(xpv)
%fpv=M\N*fcs;
fpv=linsolve(M,N*fcs); 
 
end

switch flag,
%----- S-FUNCTION BUILDING PARAMETERS AND INITIAL CONDITIONS --------
    case 0,        
sys(1)=7;%Number of continuous states
sys(2)=0;%Number of discrete states
sys(3)=7;%Number of outputs
sys(4)=7;%Number of inputs
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
%----- Process variable derivatives:
%5x5 linear system:
%vec=[dhWCSdt; dhCSdt; dhCLdt; dPdt; dTdt];
% M=[2*LCS*sqrt(hWCS*(D-hWCS))   	0       				0 						-VWCSstd*dBwdP   															-VWCSstd*dBwdT
   % -2*LCS*sqrt(hWCS*(D-hWCS))  	2*LCS*sqrt(hCS*(D-hCS)) 0 						-VWLCSstd*dBwdP-VOLCSstd*dBodP  											-VWLCSstd*dBwdT-VOLCSstd*dBodT
   % 0   							0       				2*LCL*sqrt(hCL*(D-hCL)) -VWLCLstd*dBwdP-VOLCLstd*dBodP  											-VWLCLstd*dBwdT-VOLCLstd*dBodT
   % 0 							0 						0  						VGstd*dBgdP+(VWCSstd+VWLCSstd+VWLCLstd)*dBwdP+(VOLCSstd+VOLCLstd)*dBodP   	VGstd*dBgdT+(VWCSstd+VWLCSstd+VWLCLstd)*dBwdT+(VOLCSstd+VOLCLstd)*dBodT
   % 0   							0       				0  						0 																			(VWCSstd+VWLCSstd+VWLCLstd)*rWstd*cPw+(VOLCSstd+VOLCLstd)*d60*PP.rho_H2Ostd*cPo];
% b=[Bw*dVWCSstddt; ...
   % Bw*dVWLCSstddt+Bo*dVOLCSstddt; ...
   % Bw*dVWLCLstddt+Bo*dVOLCLstddt; ...
   % -Bg*dVGstddt-Bw*(dVWCSstddt+dVWLCSstddt+dVWLCLstddt)-Bo*(dVOLCSstddt+dVOLCLstddt); ...
   % dEdt-rWstd*cPw*(T-PP.Tref)*(dVWCSstddt+dVWLCSstddt+dVWLCLstddt)-d60*PP.rho_H2Ostd*cPo*(T-PP.Tref)*(dVOLCSstddt+dVOLCLstddt)];
% vec=M\b;
% dhWCSdt=vec(1);
% dhCSdt=vec(2);
% dhCLdt=vec(3);
% dPdt=vec(4);
% dTdt=vec(5);

% %Water std. vol. content at light phase in CS rate:
% dBSWLCSdt=((1-BSWLCS)*dVWLCSstddt-BSWLCS*dVOLCSstddt)/(VOLCSstd+VWLCSstd);

% %Water std. vol. content at light phase in CL rate:
% dBSWLCLdt=((1-BSWLCL)*dVWLCLstddt-BSWLCL*dVOLCLstddt)/(VOLCLstd+VWLCLstd);

dhWCSdt		= fpv(1);
dBSWLCSdt	= fpv(2);
dhCSdt		= fpv(3);
dhCLdt		= fpv(4);
dBSWLCLdt	= fpv(5);
dPdt		= fpv(6);
dTdt		= fpv(7);

%----- State derivatives for sys output flag=1
sys(1)=dhWCSdt;
sys(2)=dBSWLCSdt;
sys(3)=dhCSdt;
sys(4)=dhCLdt;
sys(5)=dBSWLCLdt;
sys(6)=dPdt;
sys(7)=dTdt;
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
	case 10,	%dfdx
%----- Canonical states dfdx:
%M*fpv - N*fcs = 0;
%d(M*fpv)dxpv(i) - d(N*fcs)dpxv(i) = 0;
%M*dfpvdxpv(i) + dMdxpv(i)*fpv + sum[j=1:Nxcs]{dMdxcs(j)*dxcs(j)dxpv(i)*fpv} ...
% - N*dfcsdxpv(i) - dNdxpv(i)*fcs - sum[j=1:Nxcs]{dNdxcs(j)*dxcs(j)dxpv(i)*fcs} = 0;

%M*fpv - N*fcs = 0 -> M*fpv - N*fcs -> M*dxpvdt - N*dxcsdt = 0;
%dxcsdt/dxpvdt = N\M -> dxcsdxpv = N\M -> dxcsdxpv(i) = N\M(:,i);
%dcxs(j)dxpv(i) = dxcsdxpv(j,i);

%dfpvdxpv(:,i) = M\(...
%	N*dfcsdxpv(:,i) + dNdxpv(:,:,i)*fcs + sum[j=1:Nxcs]{dNdxcs(:,:,j)*dxcsdxpv(j,i)*fcs} ...
%	- dMdxpv(:,:,i)*fpv - sum[j=1:Nxcs]{dMdxcs(:,:,j)*dxcsdxpv(j,i)*fpv} ...
%	);

%xcs partial derivative in relation to xpv
dxcsdxpv=N\M;

%----- dfcsdxpv:
%Qweir=Qweir(BSWLCS,hCS,hCL,P,T);
%effWL=effWL(hWCS,Bo_in(P));
%Bo=Bo(P,T);
%Bw=Bw(P,T);
%Bg=Bg(P,T);

deffWLdhWCS=...
	-0.2675 + 0.1356*(Qf_in/(24*3600)*(1-BSW_in)*Bo_in) - 0.0747*2*(hWCS);
dfcsdxpv(:,1)=[...	%dfcsdhWCS
	0; ...
	(-Qf_in*BSW_in*deffWLdhWCS)/(24*3600); ...
	(Qf_in*BSW_in*deffWLdhWCS)/(24*3600); ...
	0; ...
	0; ...
	0; ...
	0];

dQweirdBSWLCS=-(Bw-Bo)/((1-BSWLCS)*Bo+BSWLCS*Bw)*Qweir;
dfcsdxpv(:,2)=[...	%dfcsdBSWLCS
	(Qweir-dQweirdBSWLCS*(1-BSWLCS))/(24*3600); ...
	(-Qweir-dQweirdBSWLCS*BSWLCS)/(24*3600); ...
	0; ...
	(-Qweir+dQweirdBSWLCS*(1-BSWLCS))/(24*3600); ...
	(Qweir+dQweirdBSWLCS*BSWLCS)/(24*3600); ...
	0; ...
	0];	

dQweirdhCS=(24*3600)/((1-BSWLCS)*Bo+BSWLCS*Bw)*( ...
	+ 0.5*(1+tanh((hweir-hCL)/h))*0.5*(1-tanh((hCS-hweir)/h)^2)*(1/h)* ...
		(24.88*sqrt(2*9.81)*(Lweir-0.2*(abs(hCS-hweir)))*(abs(hCS-hweir))^1.5) ...
	+ 0.5*(1+tanh((hweir-hCL)/h))*0.5*(1+tanh((hCS-hweir)/h))* ...
		(24.88*sqrt(2*9.81)*( (-0.2/sqrt((hCS-hweir)^2)*(hCS-hweir))*(abs(hCS-hweir))^1.5 + (Lweir-0.2*(abs(hCS-hweir)))*(1.5*sqrt((hCS-hweir)^2)^0.5/sqrt((hCS-hweir)^2)*(hCS-hweir)) )) ...
	+ 0.5*(1+tanh((hCL-hweir)/h))*0.5*(1-tanh((hCS-hCL)/h)^2)*(1/h)* ...
		(24.88*sqrt(2*9.81)*(2*sqrt(hCL*(D-hCL))-0.2*(abs(hCS-hCL)))*(abs(hCS-hCL))^1.5) ...
	+ 0.5*(1+tanh((hCL-hweir)/h))*0.5*(1+tanh((hCS-hCL)/h))* ...
		(24.88*sqrt(2*9.81)*( (-0.2/sqrt((hCS-hweir)^2)*(hCS-hweir))*(abs(hCS-hCL))^1.5 + (2*sqrt(hCL*(D-hCL))-0.2*(abs(hCS-hCL)))*(1.5*sqrt((hCS-hCL)^2)^0.5/sqrt((hCS-hCL)^2)*(hCS-hCL)) )) ...
	);
dfcsdxpv(:,3)=[...	%dfcsdhCS
	(-dQweirdhCS*(1-BSWLCS))/(24*3600); ...
	(-dQweirdhCS*BSWLCS)/(24*3600); ...
	0; ...
	(dQweirdhCS*(1-BSWLCS))/(24*3600); ...
	(dQweirdhCS*BSWLCS)/(24*3600); ...
	0; ...
	0];	

dQweirdhCL=(24*3600)/((1-BSWLCS)*Bo+BSWLCS*Bw)*( ...
	+ 0.5*(1-tanh((hweir-hCL)/h)^2)*(-1/h)*0.5*(1+tanh((hCS-hweir)/h))* ...
		(24.88*sqrt(2*9.81)*(Lweir-0.2*(abs(hCS-hweir)))*(abs(hCS-hweir))^1.5) ...
	+ 0.5*(1-tanh((hCL-hweir)/h)^2)*(1/h)*0.5*(1+tanh((hCS-hCL)/h))* ...
		(24.88*sqrt(2*9.81)*(2*sqrt(hCL*(D-hCL))-0.2*(abs(hCS-hCL)))*(abs(hCS-hCL))^1.5) ...
	+ 0.5*(1+tanh((hCL-hweir)/h))*0.5*(1-tanh((hCS-hCL)/h)^2)*(-1/h)* ...
		(24.88*sqrt(2*9.81)*(2*sqrt(hCL*(D-hCL))-0.2*(abs(hCS-hCL)))*(abs(hCS-hCL))^1.5) ...
	+ 0.5*(1+tanh((hCL-hweir)/h))*0.5*(1+tanh((hCS-hCL)/h))* ...
		(24.88*sqrt(2*9.81)*( (1/sqrt(hCL*(D-hCL))*(-2*hCL+D)-0.2/sqrt((hCS-hCL)^2)*(-(hCS-hCL)))*(abs(hCS-hCL))^1.5 + (2*sqrt(hCL*(D-hCL))-0.2*(abs(hCS-hCL)))*(1.5*sqrt((hCS-hCL)^2)^0.5/sqrt((hCS-hCL)^2)*(-(hCS-hCL))) )) ...		
	);
dfcsdxpv(:,4)=[...	%dfcsdhCL
	(-dQweirdhCL*(1-BSWLCS))/(24*3600); ...
	(-dQweirdhCL*BSWLCS)/(24*3600); ...
	0; ...
	(dQweirdhCL*(1-BSWLCS))/(24*3600); ...
	(dQweirdhCL*BSWLCS)/(24*3600); ...
	0; ...
	0];	

dfcsdxpv(:,5)=[...	%dfcsdBSWLCL
	0; ...
	0; ...
	0; ...
	(Qf)/(24*3600); ...
	(-Qf)/(24*3600); ...
	(Qf*Rs)/(24*3600); ...
	(Qf*d60*PP.rho_H2Ostd*cPo*(T-PP.Tref) -Qf*rWstd*cPw*(T-PP.Tref))/(24*3600)];

dQweirdP=-((1-BSWLCS)*dBodP+BSWLCS*dBwdP)/((1-BSWLCS)*Bo+BSWLCS*Bw)*Qweir;
dBo_indP=dOilFVFdP(PP,P,T_in,d60,dg);
deffWLdP=-0.3626*(Qf_in/(24*3600)*(1-BSW_in)*dBo_indP) ...
	+ 0.1356*(hWCS)*(Qf_in/(24*3600)*(1-BSW_in)*dBo_indP) ...
	+ 0.089*(2)*(Qf_in/(24*3600)*(1-BSW_in)*Bo_in)*(Qf_in/(24*3600)*(1-BSW_in)*dBo_indP);
dRs_indP=dSolGasOilRatiodP(PP,P,T_in,d60,dg); dRsdP=dSolGasOilRatiodP(PP,P,T,d60,dg);
deGres_indP=dGasSpecResEnthdP(PP,T_in,dg); deGresdP=dGasSpecResEnthdP(PP,T,dg);
dfcsdxpv(:,6)=[...	%dfcsdP
	(-Qf_in*BSW_in*deffWLdP-dQweirdP*(1-BSWLCS))/(24*3600); ...
	(Qf_in*BSW_in*deffWLdP-dQweirdP*BSWLCS)/(24*3600); ...
	0; ...
	(dQweirdP*(1-BSWLCS))/(24*3600); ...
	(dQweirdP*BSWLCS)/(24*3600); ...
	(Qf_in*(1-BSW_in)*dRs_indP-Qf*(1-BSWLCL)*dRsdP)/(24*3600); ...
	(Qg_in*dg*PP.rho_AIRstd*deGres_indP-Qg*dg*PP.rho_AIRstd*deGresdP)/(24*3600)];

dQweirdT=-((1-BSWLCS)*dBodT+BSWLCS*dBwdT)/((1-BSWLCS)*Bo+BSWLCS*Bw)*Qweir;
dRsdT=dSolGasOilRatiodT(PP,P,T,d60,dg);
deGresdT=dGasSpecResEnthdT(PP,P,T,dg);
dfcsdxpv(:,7)=[...	%dfcsdT
	(-dQweirdT*(1-BSWLCS))/(24*3600); ...
	(-dQweirdT*BSWLCS)/(24*3600); ...
	0; ...
	(dQweirdT*(1-BSWLCS))/(24*3600); ...
	(dQweirdT*BSWLCS)/(24*3600); ...
	(-Qf*(1-BSWLCL)*dRsdT)/(24*3600); ...
	(-Qf*(1-BSWLCL)*d60*PP.rho_H2Ostd*cPo-Qf*BSWLCL*rWstd*cPw-Qg*dg*PP.rho_AIRstd*(deGresdT+cPg)-Qw*rWstd*cPw)/(24*3600)];

%----- dMdxcs:
dMdxcs(:,:,1)=[...	%dMdVOLCSstd
	0   	0		0       0 		0		0   	0
	0		eye(1)	0		0		0		0		0
	0  		0		0 		0 		0		-dBodP  -dBodT
	0   	0		0       0 		0		0  		0
	0		0		0		0		0		0		0
	0 		0		0 		0  		0		dBodP   dBodT
	0   	0		0       0  		0		0 		d60*PP.rho_H2Ostd*cPo];
	
dMdxcs(:,:,2)=[...	%dMdVWLCSstd
	0   	0		0       0 		0		0   	0
	0		eye(1)	0		0		0		0		0
	0  		0		0 		0 		0		-dBwdP  -dBwdT
	0   	0		0       0 		0		0  		0
	0		0		0		0		0		0		0
	0 		0		0 		0  		0		dBwdP   dBwdT
	0   	0		0       0  		0		0 		rWstd*cPw];

dMdxcs(:,:,3)=[...	%dMdVWCSstd
	0   	0		0       0 		0		-dBwdP   -dBwdT
	0		0		0		0		0		0		0
	0  		0		0 		0 		0		0  		0
	0   	0		0       0 		0		0  		0
	0		0		0		0		0		0		0
	0 		0		0 		0  		0		dBwdP   dBwdT
	0   	0		0       0  		0		0 		rWstd*cPw];

dMdxcs(:,:,4)=[...	%dMdVOLCLstd
	0   	0		0       0 		0		0   	0
	0		0		0		0		0		0		0
	0  		0		0 		0 		0		0  		0
	0   	0		0       0 		0		-dBodP  -dBodT
	0		0		0		0		eye(1)	0		0
	0 		0		0 		0  		0		dBodP   dBodT
	0   	0		0       0  		0		0 		d60*PP.rho_H2Ostd*cPo];

dMdxcs(:,:,5)=[...	%dMdVWLCLstd
	0   	0		0       0 		0		0   	0
	0		0		0		0		0		0		0
	0  		0		0 		0 		0		0  		0
	0   	0		0       0 		0		-dBwdP  -dBwdT
	0		0		0		0		eye(1)	0		0
	0 		0		0 		0  		0		dBwdP   dBwdT
	0   	0		0       0  		0		0 		rWstd*cPw];

dMdxcs(:,:,6)=[...	%dMdVGstd
	0   	0		0       0 		0		0   	0
	0		0		0		0		0		0		0
	0  		0		0 		0 		0		0  		0
	0   	0		0       0 		0		0  		0
	0		0		0		0		0		0		0
	0 		0		0 		0  		0		dBgdP   dBgdT
	0   	0		0       0  		0		0 		0];

dMdxcs(:,:,7)=[...	%dMdE
	0   	0		0       0 		0		0   	0
	0		0		0		0		0		0		0
	0  		0		0 		0 		0		0  		0
	0   	0		0       0 		0		0  		0
	0		0		0		0		0		0		0
	0 		0		0 		0  		0		0   	0
	0   	0		0       0  		0		0 		0];

%----- dMdxpv:
d2BodP2=d2OilFVFdP2(PP,P,T,d60,dg); d2BodPdT=d2OilFVFdPdT(PP,P,T,d60,dg);
d2BodTdP=d2OilFVFdTdP(PP,P,T,d60,dg); d2BodT2=d2OilFVFdT2(PP,P,T,d60,dg);

d2BwdP2=d2WaterFVFdP2(PP,T); d2BwdPdT=d2WaterFVFdPdT(PP,P,T);
d2BwdTdP=d2WaterFVFdTdP(PP,P,T); d2BwdT2=d2WaterFVFdT2(PP,P,T);

d2BgdP2=d2GasFVFdP2(PP,P,T,dg); d2BgdPdT=d2GasFVFdPdT(PP,P,T,dg);
d2BgdTdP=d2GasFVFdTdP(PP,P,T,dg); d2BgdT2=d2GasFVFdT2(PP,P,T,dg);


dMdxpv(:,:,1)=[...	%dMdhWcs
	2*LCS*(1/2)/sqrt(hWCS*(D-hWCS))*(-2*hWCS+D)   	0						0       									0 											0						0   																				0
	0												0						0											0											0						0																					0
	-2*LCS*(1/2)/sqrt(hWCS*(D-hWCS))*(-2*hWCS+D)  	0						0 											0 											0						0  																					0
	0   											0						0       									0 											0						0  																					0
	0												0						0											0											0						0																					0
	0 												0						0 											0  											0						0   																				0
	0   											0						0       									0  											0						0 																					0];

dMdxpv(:,:,2)=[...	%dMdBSWLCS							
	0   											0						0       									0 											0						0   																				0
	0												0						0											0											0						0																					0
	0  												0						0 											0 											0						0  																					0
	0   											0						0       									0 											0						0  																					0
	0												0						0											0											0						0																					0
	0 												0						0 											0  											0						0   																				0
	0   											0						0       									0  											0						0 																					0];

dMdxpv(:,:,3)=[...	%dMdhCS							
	0   											0						0       									0 											0						0   																				0
	0												0						0											0											0						0																					0
	0  												0						2*LCS*(1/2)/sqrt(hCS*(D-hCS))*(-2*hCS+D) 	0 											0						0  																					0
	0   											0						0       									0 											0						0  																					0
	0												0						0											0											0						0																					0
	0 												0						0 											0  											0						0   																				0
	0   											0						0       									0  											0						0 																					0];
		
dMdxpv(:,:,4)=[...	%dMdhCL		
	0   											0						0       									0 											0						0   																				0
	0												0						0											0											0						0																					0
	0  												0						0 											0 											0						0  																					0
	0   											0						0       									2*LCL*(1/2)/sqrt(hCL*(D-hCL))*(-2*hCL+D) 	0						0  																					0
	0												0						0											0											0						0																					0
	0 												0						0 											0  											0						0   																				0
	0   											0						0       									0  											0						0 																					0];
		
dMdxpv(:,:,5)=[...	%dMdBSWLCL		
	0   											0						0       									0 											0						0   																				0
	0												0						0											0											0						0																					0
	0  												0						0 											0 											0						0  																					0
	0   											0						0       									0 											0						0  																					0
	0												0						0											0											0						0																					0
	0 												0						0 											0  											0						0   																				0
	0   											0						0       									0  											0						0 																					0];

dMdxpv(:,:,6)=[...	%dMdP
	0   											0						0       									0 											0						-VWCSstd*d2BwdP2   																	-VWCSstd*d2BwdTdP
	0												0						0											0											0						0																					0
	0  												0						0 											0 											0						-VWLCSstd*d2BwdP2-VOLCSstd*d2BodP2  												-VWLCSstd*d2BwdTdP-VOLCSstd*d2BodTdP
	0   											0						0       									0 											0						-VWLCLstd*d2BwdP2-VOLCLstd*d2BodP2  												-VWLCLstd*d2BwdTdP-VOLCLstd*d2BodTdP
	0												0						0											0											0						0																					0
	0 												0						0 											0  											0						VGstd*d2BgdP2+(VWCSstd+VWLCSstd+VWLCLstd)*d2BwdP2+(VOLCSstd+VOLCLstd)*d2BodP2   	VGstd*d2BgdTdP+(VWCSstd+VWLCSstd+VWLCLstd)*d2BwdTdP+(VOLCSstd+VOLCLstd)*d2BodTdP
	0   											0						0       									0  											0						0 																					0];

dMdxpv(:,:,7)=[...	%dMdT
	0   											0						0       									0 											0						-VWCSstd*d2BwdPdT   																-VWCSstd*d2BwdT2
	0												0						0											0											0						0																					0
	0  												0						0 											0 											0						-VWLCSstd*d2BwdPdT-VOLCSstd*d2BodPdT  												-VWLCSstd*d2BwdT2-VOLCSstd*d2BodT2
	0   											0						0       									0 											0						-VWLCLstd*d2BwdPdT-VOLCLstd*d2BodPdT  												-VWLCLstd*d2BwdT2-VOLCLstd*d2BodT2
	0												0						0											0											0						0																					0
	0 												0						0 											0  											0						VGstd*d2BgdPdT+(VWCSstd+VWLCSstd+VWLCLstd)*d2BwdPdT+(VOLCSstd+VOLCLstd)*d2BodPdT   	VGstd*d2BgdT2+(VWCSstd+VWLCSstd+VWLCLstd)*d2BwdT2+(VOLCSstd+VOLCLstd)*d2BodT2
	0   											0						0       									0  											0						0 																					0];

%----- dNdxcs:
dNdxcs(:,:,1)=[...	%dNdVOLCSstd
	0   	0		0       0 		0		0   	0
	0		0		0		0		0		0		0
	0  		0		0 		0 		0		0  		0
	0   	0		0       0 		0		0  		0
	0		0		0		0		0		0		0
	0 		0		0 		0  		0		0   	0
	0   	0		0       0  		0		0 		0];
	
dNdxcs(:,:,2)=[...	%dNdVWLCSstd
	0   	0		0       0 		0		0   	0
	0		0		0		0		0		0		0
	0  		0		0 		0 		0		0  		0
	0   	0		0       0 		0		0  		0
	0		0		0		0		0		0		0
	0 		0		0 		0  		0		0   	0
	0   	0		0       0  		0		0 		0];

dNdxcs(:,:,3)=[...	%dNdVWCSstd
	0   	0		0       0 		0		0   	0
	0		0		0		0		0		0		0
	0  		0		0 		0 		0		0  		0
	0   	0		0       0 		0		0  		0
	0		0		0		0		0		0		0
	0 		0		0 		0  		0		0   	0
	0   	0		0       0  		0		0 		0];

dNdxcs(:,:,4)=[...	%dNdVOLCLstd
	0   	0		0       0 		0		0   	0
	0		0		0		0		0		0		0
	0  		0		0 		0 		0		0  		0
	0   	0		0       0 		0		0  		0
	0		0		0		0		0		0		0
	0 		0		0 		0  		0		0   	0
	0   	0		0       0  		0		0 		0];

dNdxcs(:,:,5)=[...	%dNdVWLCLstd
	0   	0		0       0 		0		0   	0
	0		0		0		0		0		0		0
	0  		0		0 		0 		0		0  		0
	0   	0		0       0 		0		0  		0
	0		0		0		0		0		0		0
	0 		0		0 		0  		0		0   	0
	0   	0		0       0  		0		0 		0];

dNdxcs(:,:,6)=[...	%dNdVGstd
	0   	0		0       0 		0		0   	0
	0		0		0		0		0		0		0
	0  		0		0 		0 		0		0  		0
	0   	0		0       0 		0		0  		0
	0		0		0		0		0		0		0
	0 		0		0 		0  		0		0   	0
	0   	0		0       0  		0		0 		0];

dNdxcs(:,:,7)=[...	%dNdE
	0   	0		0       0 		0		0   	0
	0		0		0		0		0		0		0
	0  		0		0 		0 		0		0  		0
	0   	0		0       0 		0		0  		0
	0		0		0		0		0		0		0
	0 		0		0 		0  		0		0   	0
	0   	0		0       0  		0		0 		0];

%----- dNdxpv:
dNdxpv(:,:,1)=[...	%dNdhWCS
	0   						0				0       		0 								0				0   	0
	0							0				0				0								0				0		0
	0  							0				0 				0 								0				0  		0
	0   						0				0       		0 								0				0  		0
	0							0				0				0								0				0		0
	0 							0				0 				0  								0				0   	0
	0   						0				0       		0  								0				0 		0];


dNdxpv(:,:,2)=[...	%dNdVBSWLCS												
	0   						0				0       		0 								0				0   	0
	-eye(1)						-eye(1)			0				0								0				0		0
	0  							0				0 				0 								0				0  		0
	0   						0				0       		0 								0				0  		0
	0							0				0				0								0				0		0
	0 							0				0 				0  								0				0   	0
	0   						0				0       		0  								0				0 		0];

dNdxpv(:,:,3)=[...	%dNdhCS												
	0   						0				0       		0 								0				0   	0
	0							0				0				0								0				0		0
	0  							0				0 				0 								0				0  		0
	0   						0				0       		0 								0				0  		0
	0							0				0				0								0				0		0
	0 							0				0 				0  								0				0   	0
	0   						0				0       		0  								0				0 		0];

dNdxpv(:,:,4)=[...	%dNdVhCL												
	0   						0				0       		0 								0				0   	0
	0							0				0				0								0				0		0
	0  							0				0 				0 								0				0  		0
	0   						0				0       		0 								0				0  		0
	0							0				0				0								0				0		0
	0 							0				0 				0  								0				0   	0
	0   						0				0       		0  								0				0 		0];

dNdxpv(:,:,5)=[...	%dNdBSWLCL												
	0   						0				0       		0 								0				0   	0
	0							0				0				0								0				0		0
	0  							0				0 				0 								0				0  		0
	0   						0				0       		0 								0				0  		0
	0							0				0				-eye(1)							-eye(1)			0		0
	0 							0				0 				0  								0				0   	0
	0   						0				0       		0  								0				0 		0];

dNdxpv(:,:,6)=[...	%dNdP										
	0   						0				dBwdP			0 								0				0   	0
	0							0				0				0								0				0		0
	dBodP  						dBwdP			0 				0 								0				0  		0
	0   						0				0       		dBodP 							dBwdP			0  		0
	0							0				0				0								0				0		0
	-dBodP 						-dBwdP			-dBwdP 			-dBodP  						-dBwdP			-dBgdP	0
	0   						0				0       		0  								0				0 		0];

dNdxpv(:,:,7)=[...	%dNdT
	0   						0				dBwdT			0 								0				0   	0
	0							0				0				0								0				0		0
	dBodT  						dBwdT			0 				0 								0				0  		0
	0   						0				0       		dBodT 							dBwdT			0  		0
	0							0				0				0								0				0		0
	-dBodT 						-dBwdT			-dBwdT 			-dBodT  						-dBwdT			-dBgdT	0
	-d60*PP.rho_H2Ostd*cPo   	-rWstd*cPw		-rWstd*cPw      -d60*PP.rho_H2Ostd*cPo  		-rWstd*cPw		0 		0];

%----- dfdx for sys output flag=10
for i=1:7
	dfdx(:,i) = M\(N*dfcsdxpv(:,i) + dNdxpv(:,:,i)*fcs - dMdxpv(:,:,i)*fpv);
	for j=1:7
		dfdx(:,i) = dfdx(:,i) + M\(...
			dNdxcs(:,:,j)*dxcsdxpv(j,i)*fcs - dMdxcs(:,:,j)*dxcsdxpv(j,i)*fpv);
	end
end

sys=dfdx;
%--------------------------------------------------------------------	
	case 11,	%dfdu
%----- Canonical states dfdu:
%M=dhdxpv(xcs,xpv,par); N=-dhdxcs(xcs,xpv,par);
%M*fpv - N*fcs = 0;
%d(M*fpv)du(i) - d(N*fcs)du(i) = 0;
%M*dfpvdu - N*dfcsdu = 0 -> dfpvdu=M\N*dfcsdu

%----- dfcsdu:	
%effWL=effWL(Qf_in,BSW_in,Bo_in(T_in));

deffWLdQf_in=-0.3626*(1/(24*3600)*(1-BSW_in)*Bo_in) ...
	+ 0.1356*(hWCS)*(1/(24*3600)*(1-BSW_in)*Bo_in) ...
	+ 0.089*(2)*(Qf_in/(24*3600)*(1-BSW_in)*Bo_in)*(1/(24*3600)*(1-BSW_in)*Bo_in);
dfcsdu(:,1)=[...	%dfcsdQf_in
	((1-BSW_in))/(24*3600); ...
	(BSW_in*(1-effWL) - Qf_in*BSW_in*deffWLdQf_in)/(24*3600); ...
	(BSW_in*effWL + Qf_in*BSW_in*deffWLdQf_in)/(24*3600); ...
	0; ...
	0; ...
	((1-BSW_in)*Rs_in)/(24*3600); ...
	((1-BSW_in)*d60*PP.rho_H2Ostd*cPo*(T_in-PP.Tref) + BSW_in*rWstd*cPw*(T_in-PP.Tref))/(24*3600)];

dfcsdu(:,2)=[...	%dfcsdQg_in
	0; ...
	0; ...
	0; ...
	0; ...
	0; ...
	1/(24*3600); ...
	(dg*PP.rho_AIRstd*(eGres_in+DeOG+cPg*(T_in-PP.Tref)))/(24*3600)];

deffWLdBSW_in=-0.3626*(Qf_in/(24*3600)*(-1)*Bo_in) ...
	+ 0.1356*(hWCS)*(Qf_in/(24*3600)*(-1)*Bo_in) ...
	+ 0.089*(2)*(Qf_in/(24*3600)*(1-BSW_in)*Bo_in)*(Qf_in/(24*3600)*(-1)*Bo_in);
dfcsdu(:,3)=[...	%dfcsdBSW_in
	(-Qf_in)/(24*3600); ...
	(Qf_in*(1-effWL) - Qf_in*BSW_in*deffWLdBSW_in)/(24*3600); ...
	(Qf_in*effWL + Qf_in*BSW_in*deffWLdBSW_in)/(24*3600); ...
	0; ...
	0; ...
	(-Qf_in*Rs_in)/(24*3600); ...
	(-Qf_in*d60*PP.rho_H2Ostd*cPo*(T_in-PP.Tref) + Qf_in*rWstd*cPw*(T_in-PP.Tref))/(24*3600)];	

dBo_indT_in=dOilFVFdT(PP,P,T_in,d60,dg);
deffWLdT_in=-0.3626*(Qf_in/(24*3600)*(1-BSW_in)*dBo_indT_in) ...
	+ 0.1356*(hWCS)*(Qf_in/(24*3600)*(1-BSW_in)*dBo_indT_in) ...
	+ 0.089*(2)*(Qf_in/(24*3600)*(1-BSW_in)*Bo_in)*(Qf_in/(24*3600)*(1-BSW_in)*dBo_indT_in);
dRs_indT_in=dSolGasOilRatiodT(PP,P,T_in,d60,dg);	
deGres_indT_in=dGasSpecResEnthdT(PP,P,T_in,dg);
dfcsdu(:,4)=[...	%dfcsdT_in
	0; ...
	(-Qf_in*BSW_in*deffWLdT_in)/(24*3600); ...
	(Qf_in*BSW_in*deffWLdT_in)/(24*3600); ...
	0; ...
	0; ...
	(Qf_in*(1-BSW_in)*dRs_indT_in)/(24*3600); ...
	(Qf_in*(1-BSW_in)*d60*PP.rho_H2Ostd*cPo + Qf_in*BSW_in*rWstd*cPw + Qg_in*dg*PP.rho_AIRstd*(deGres_indT_in+cPg))/(24*3600)];

dfcsdu(:,5)=[...	%dfcsdQf
	0; ...
	0; ...
	0; ...
	(-(1-BSWLCL))/(24*3600); ...
	(-BSWLCL)/(24*3600); ...
	(-(1-BSWLCL)*Rs)/(24*3600); ...
	(-(1-BSWLCL)*d60*PP.rho_H2Ostd*cPo*(T-PP.Tref) - BSWLCL*rWstd*cPw*(T-PP.Tref))/(24*3600)];

dfcsdu(:,6)=[...	%dfcsdQw
	0; ...
	0; ...
	(-1)/(24*3600); ...
	0; ...
	0; ...
	0; ...
	(-rWstd*cPw*(T-PP.Tref))/(24*3600)];

dfcsdu(:,7)=[...	%dfcsdQg
	0; ...
	0; ...
	0; ...
	0; ...
	0; ...
	(-1)/(24*3600); ...
	(-dg*PP.rho_AIRstd*(eGres+DeOG+cPg*(T-PP.Tref)))/(24*3600)];
	
%----- dfdu for sys output flag=11	
dfdu=M\N*dfcsdu;
sys=dfdu;
%--------------------------------------------------------------------	
	case 12,	%dhdx
dhdx=eye(7,7);
sys=dhdx;
%--------------------------------------------------------------------	
	case 13,	%dhdu
dhdu=zeros(7,7);
sys=dhdu;
%--------------------------------------------------------------------	
    otherwise   % Unexpected flags
DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
end
%--------------------------------------------------------------------    
end