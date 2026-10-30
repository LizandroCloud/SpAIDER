function [sys,x0,str,ts,simStateCompliance]=scrubber_gas(t,x,u,flag,cinit,par)
%scrubber_gas - Gas-only scrubber vessel non-linear model S-Function
%
%   This model describes the dynamic behavior of an upstream oil &
%   gas process scrubber vessel assuming gas phase only, using the 
%	black-oil approach.
%
%   [sys,x0,str,ts,simStateCompliance]=scrubber_gas(t,x,u,flag,cinit,par)
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
%   par{1}: d60     Oil specific gravity              		 []
%   par{2}: dg      Gas specific gravity                     []
%   par{3}: C       Vessel length                            [m]
%   par{4}: D       Vessel diameter                          [m]
%
%   Canonical states:
%   xcs(1): VGstd     Gas std. vol. hold-up                  [m^3]
%   xcs(2): E         Energy hold-up                         [kJ]
%
%   Process variable states:
%   xpv(1): P         Pressure                               [kPa]
%   xpv(2): T         Temperature                            [K]
%
%   Inputs:
%   u(1): Qg_in     Inlet gas std. vol. flow rate            [Sm^3/d]
%   u(2): T_in      Inlet temperature                        [K]
%   u(3): Qg        Gas std. vol. flow rate                  [Sm^3/d]
%
%   Outputs:
%   y(1): P         Pressure                                 [kPa]
%   y(2): T         Temperature                              [K]
%
%   Initial conditions set:
%   cinit = xpv

%   2023. Coded by Diego Trica (diegotrica@gmail.com)

%Blackoil object class property package. On the main script the
%user must declare PP as global and PP=blackoil.
    global PP

%----- Parameters:
d60     = par{1};   %Oil specific gravity []
dg      = par{2};   %Gas specific gravity []
C       = par{3};   %Vessel length [m]
D       = par{4};   %Vessel diameter [m]

%Vessel total volume [m^3]
Vt=0.25*C*pi*D^2;

%Gas specific heat capacity [kJ/kg K]
cPg=GasSpecHeatCap(PP,dg);

%Gas liberation specific enthapy [kJ/kg]
DeOG=GasLatHeat(PP,d60,dg);
   
%----- States:
P 		= x(1);
T		= x(2);

%----- Inputs:
Qg_in	= u(1);
T_in	= u(2);
Qg      = u(3);   

if flag==1 || flag==10 || flag==11
%----- Internal calculations:
%Gas formation volume factor [m^3/Sm^3]
Bg=GasFVF(PP,P,T,dg);

%Inlet gas specific residual enthalpy [kJ/kg]
eGres_in=GasSpecResEnth(PP,P,T_in,dg);

%Outlet gas specific residual enthalpy [kJ/kg]
eGres=GasSpecResEnth(PP,P,T,dg);

%Gas FVF derivative to P and T
dBgdP=dGasFVFdP(PP,P,T,dg);
dBgdT=dGasFVFdT(PP,P,T,dg);

%gas specific residual enthalpy derivative to P and T
deGresdP=dGasSpecResEnthdP(PP,T,dg);
deGresdT=dGasSpecResEnthdT(PP,P,T,dg);

%----- Canonical states:
%Gas phase std. vol. hold-up [Sm^3]
VGstd=Vt/Bg;

%----- Canonical state derivatives:
%Gas std. vol. balance
dVGstddt=(Qg_in-Qg)/(24*3600);

%Energy balance
dEdt=(24*3600)\(...
 +Qg_in*dg*PP.rho_AIRstd*(eGres_in+DeOG+cPg*(T_in-PP.Tref))...
 -Qg*dg*PP.rho_AIRstd*(eGres+DeOG+cPg*(T-PP.Tref)));

%----- Process variable derivatives:
%h(xcs,xpv,par)=0;	%algebraic relation between xpv and xcs
%dhdt = 0 -> dhdxpv*dxpvdt + dhdxcs*dxcsdt = 0;
%dxcsdt=fcs(xcs,xpv,u,par); dxpvdt=fpv(xcs,xpv,u,par);
%dhdxpv*fpv + dhdxcs*fcs = 0;
%M=dhdxpv(xcs,xpv,par); N=-dhdxcs(xcs,xpv,par);
%M*fpv - N*fcs = 0 -> fpv = M\N*fcs;

%dhdxpv: Algebraic relations between xcs and xpv partial derivative in relation to xpv
M=[VGstd*dBgdP   									VGstd*dBgdT
   VGstd*dg*PP.rho_AIRstd*deGresdP 					VGstd*dg*PP.rho_AIRstd*(deGresdT+cPg)];

%dhdxcs: Algebraic relations between xcs and xpv partial derivative in relation to xcs
N=[-Bg												0
   -dg*PP.rho_AIRstd*(eGres+DeOG+cPg*(T-PP.Tref))	eye(1)]; 

%diff(xcs)
fcs=[...
	dVGstddt; ...
	dEdt];

%diff(xpv)
%fpv=M\N*fcs;
fpv=linsolve(M,N*fcs);

end

switch flag,
%----- S-FUNCTION BUILDING PARAMETERS AND INITIAL CONDITIONS --------
    case 0,        
sys(1)=2;%Number of continuous states
sys(2)=0;%Number of discrete states
sys(3)=2;%Number of outputs
sys(4)=3;%Number of inputs
sys(5)=0;%Reserved for root finding. Must be zero.
sys(6)=1;%Direct Feedthrough S-function parameter
sys(7)=1;%Number of sample times (at least one sample time is needed)

%The following section computes consistent initial condtions for the
%states for a fixed set of output and input values given by css.

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
%----- Process variable state derivatives:
%2x2 linear system:
%vec=[dPdt; dTdt];
% M=[VGstd*dBgdP   						VGstd*dBgdT
   % VGstd*dg*PP.rho_AIRstd*deGresdP 		VGstd*dg*PP.rho_AIRstd*(deGresdT+cPg)];
% b=[-Bg*dVGstddt; ...
   % dEdt-dg*PP.rho_AIRstd*(eGres+DeOG+cPg*(T-PP.Tref))*dVGstddt];
% vec=M\b;
% dPdt=vec(1);
% dTdt=vec(2);

dPdt		= fpv(1);
dTdt		= fpv(2);

%----- State derivatives for sys output flag=1
sys(1)=dPdt;
sys(2)=dTdt;
sys=sys(:);
%--------------------------------------------------------------------
    case 2,
sys=[];
%----- OUPUTS -------------------------------------------------------
    case 3,
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
%Bg=Bg(P,T);

deGres_indP=dGasSpecResEnthdP(PP,T_in,dg); deGresdP=dGasSpecResEnthdP(PP,T,dg);
dfcsdxpv(:,1)=[...	%dfcsdP
	0; ...
	(Qg_in*dg*PP.rho_AIRstd*deGres_indP - Qg*dg*PP.rho_AIRstd*deGresdP)/(24*3600)];
	
deGresdT=dGasSpecResEnthdT(PP,P,T,dg);
dfcsdxpv(:,2)=[...	%dfcsdT
	0; ...
	(-Qg*dg*PP.rho_AIRstd*(deGresdT+cPg))/(24*3600)];	

%----- dMdxcs:
dMdxcs(:,:,1)=[...	%dMdVGstd
	dBgdP   					dBgdT
	dg*PP.rho_AIRstd*deGresdP 	dg*PP.rho_AIRstd*(deGresdT+cPg)];

dMdxcs(:,:,2)=[...	%dMdE	
	0   						0	
	0 							0];

%----- dMdxpv:
d2BgdP2=d2GasFVFdP2(PP,P,T,dg); d2BgdPdT=d2GasFVFdPdT(PP,P,T,dg);
d2BgdTdP=d2GasFVFdTdP(PP,P,T,dg); d2BgdT2=d2GasFVFdT2(PP,P,T,dg);

dMdxpv(:,:,1)=[...	%dMdP
	VGstd*d2BgdP2	VGstd*d2BgdTdP
	0 				0];

dMdxpv(:,:,2)=[...	%dMdT
	VGstd*d2BgdPdT	VGstd*d2BgdT2
	0 				0];	

%----- dNdxcs:
dNdxcs(:,:,1)=[...	%dNdVGstd
	0   	0
	0 		0];

dNdxcs(:,:,2)=[...	%dNdE
	0   	0
	0 		0];

%----- dNdxpv:
dNdxpv(:,:,1)=[...	%dNdP
	-dBgdP								0
	-dg*PP.rho_AIRstd*(deGresdP)		0];

dNdxpv(:,:,2)=[...	%dNdT
	-dBgdT								0
	-dg*PP.rho_AIRstd*(deGresdT+cPg) 	0];

%----- dfdx for sys output flag=10
for i=1:2
	dfdx(:,i) = M\(N*dfcsdxpv(:,i) + dNdxpv(:,:,i)*fcs - dMdxpv(:,:,i)*fpv);
	for j=1:2
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
dfcsdu(:,1)=[...	%dfcsdQg_in
	1/(24*3600); ...
	(dg*PP.rho_AIRstd*(eGres_in+DeOG+cPg*(T_in-PP.Tref)))/(24*3600)];

deGres_indT_in=dGasSpecResEnthdT(PP,P,T_in,dg);
dfcsdu(:,2)=[...	%dfcsdT_in
	0; ...
	(Qg_in*dg*PP.rho_AIRstd*(deGres_indT_in+cPg))/(24*3600)];

dfcsdu(:,3)=[...	%dfcsdQg
	(-1)/(24*3600); ...
	(-dg*PP.rho_AIRstd*(eGres+DeOG+cPg*(T-PP.Tref)))/(24*3600)];

%----- dfdu for sys output flag=11	
dfdu=M\N*dfcsdu;
sys=dfdu;
%--------------------------------------------------------------------
    case 12,	%dhdx
dhdx=eye(2,2);
sys=dhdx;	
%--------------------------------------------------------------------
    case 13,	%dhdu
dhdu=zeros(2,3);
sys=dhdu;
%--------------------------------------------------------------------
    otherwise   % Unexpected flags
DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
end
%--------------------------------------------------------------------    
end