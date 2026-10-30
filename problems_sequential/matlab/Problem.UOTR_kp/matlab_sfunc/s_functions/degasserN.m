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

%----- Parameters:
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

if flag==1 || flag==10 || flag==11
%----- Internal calculations:
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
H=[Bo  Bw
   BSW (BSW-1)];
b=[C*Aliq; 0];
vec=H\b;
VOLstd=vec(1);
VWLstd=vec(2);

%Gas phase std. vol. hold-up [Sm^3]
VGstd=(Vt-C*Aliq)/Bg;

%----- Canonical state derivatives:
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
 
%----- Process variable derivatives:
%h(xcs,xpv,par)=0;	%algebraic relation between xpv and xcs
%dhdt = 0 -> dhdxpv*dxpvdt + dhdxcs*dxcsdt = 0;
%dxcsdt=fcs(xcs,xpv,u,par); dxpvdt=fpv(xcs,xpv,u,par);
%dhdxpv*fpv + dhdxcs*fcs = 0;
%M=dhdxpv(xcs,xpv,par); N=-dhdxcs(xcs,xpv,par);
%M*fpv - N*fcs = 0 -> fpv = M\N*fcs;

%dhdxpv: Algebraic relations between xcs and xpv partial derivative in relation to xpv
M=[2*C*sqrt(hL*(D-hL)) 						0						-VWLstd*dBwdP-VOLstd*dBodP  				-VWLstd*dBwdT-VOLstd*dBodT
   0										(VOLstd+VWLstd)			0											0
   0  										0						VGstd*dBgdP+VWLstd*dBwdP+VOLstd*dBodP   	VGstd*dBgdT+VWLstd*dBwdT+VOLstd*dBodT
   0  										0						0 											VWLstd*rWstd*cPw+VOLstd*d60*PP.rho_H2Ostd*cPo];

%dhdxcs: Algebraic relations between xcs and xpv partial derivative in relation to xcs
N=[Bo										Bw						0											0
   -BSW										(1-BSW)					0											0
   -Bo										-Bw						-Bg											0
   -d60*PP.rho_H2Ostd*cPo*(T-PP.Tref)		-rWstd*cPw*(T-PP.Tref)	0											eye(1)]; 

%diff(xcs)
fcs=[...
	dVOLstddt; ...
	dVWLstddt; ...
	dVGstddt; ...
	dEdt];

%diff(xpv)
%fpv=M\N*fcs;
fpv=linsolve(M,N*fcs);
 
end

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
%----- Process variable derivatives:
%3x3 linear system:
%vec=[dhLdt; dPdt; dTdt];
% M=[2*C*sqrt(hL*(D-hL)) 	-VWLstd*dBwdP-VOLstd*dBodP  				-VWLstd*dBwdT-VOLstd*dBodT
   % 0  					VGstd*dBgdP+VWLstd*dBwdP+VOLstd*dBodP   	VGstd*dBgdT+VWLstd*dBwdT+VOLstd*dBodT
   % 0  					0 											VWLstd*rWstd*cPw+VOLstd*d60*PP.rho_H2Ostd*cPo];
% b=[Bw*dVWLstddt+Bo*dVOLstddt; ...
   % -Bg*dVGstddt-Bw*dVWLstddt-Bo*dVOLstddt; ...
   % dEdt-rWstd*cPw*(T-PP.Tref)*dVWLstddt-d60*PP.rho_H2Ostd*cPo*(T-PP.Tref)*dVOLstddt];
% vec=M\b;
% dhLdt=vec(1);
% dPdt=vec(2);
% dTdt=vec(3);

% %Water std. vol. content at light phase in CS rate:
% dBSWdt=((1-BSW)*dVWLstddt-BSW*dVOLstddt)/(VOLstd+VWLstd);

dhLdt		= fpv(1);
dBSWdt		= fpv(2);
dPdt		= fpv(3);
dTdt		= fpv(4);

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
%Bo=Bo(P,T);
%Bw=Bw(P,T);
%Bg=Bg(P,T);

dfcsdxpv(:,1)=[...	%dfcsdhL
	0; ...
	0; ...
	0; ...
	0];

dfcsdxpv(:,2)=[...	%dfcsdBSW
	(Qf)/(24*3600); ...
	(-Qf)/(24*3600); ...
	(Qf*Rs)/(24*3600); ...
	(Qf*d60*PP.rho_H2Ostd*cPo*(T-PP.Tref) - Qf*rWstd*cPw*(T-PP.Tref))/(24*3600)];

dBo_indP=dOilFVFdP(PP,P,T_in,d60,dg);
dRs_indP=dSolGasOilRatiodP(PP,P,T_in,d60,dg); dRsdP=dSolGasOilRatiodP(PP,P,T,d60,dg);
deGres_indP=dGasSpecResEnthdP(PP,T_in,dg); deGresdP=dGasSpecResEnthdP(PP,T,dg);
dfcsdxpv(:,3)=[...	%dfcsdP
	0; ...
	0; ...
	(Qf_in*(1-BSW_in)*dRs_indP - Qf*(1-BSW)*dRsdP)/(24*3600); ...
	(Qg_in*dg*PP.rho_AIRstd*deGres_indP - Qg*dg*PP.rho_AIRstd*deGresdP)/(24*3600)];
	
dRsdT=dSolGasOilRatiodT(PP,P,T,d60,dg);
deGresdT=dGasSpecResEnthdT(PP,P,T,dg);
dfcsdxpv(:,4)=[...	%dfcsdT
	0; ...
	0; ...
	(-Qf*(1-BSW)*dRsdT)/(24*3600); ...
	(-Qf*(1-BSW)*d60*PP.rho_H2Ostd*cPo - Qf*BSW*rWstd*cPw - Qg*dg*PP.rho_AIRstd*(deGresdT+cPg))/(24*3600)];	

%----- dMdxcs:
dMdxcs(:,:,1)=[...	%dMdVOLstd
	0 		0		-dBodP  -dBodT
	0		eye(1)	0		0
	0  		0		dBodP   dBodT
	0  		0		0 		d60*PP.rho_H2Ostd*cPo];

dMdxcs(:,:,2)=[...	%dMdVWLstd
	0 		0		-dBwdP  -dBwdT
	0		eye(1)	0		0
	0  		0		dBwdP   dBwdT
	0  		0		0 		rWstd*cPw];
	
dMdxcs(:,:,3)=[...	%dMdVGstd
	0 		0		0  		0
	0		0		0		0
	0  		0		dBgdP   dBgdT
	0  		0		0 		0];

dMdxcs(:,:,4)=[...	%dMdE	
	0 		0		0  		0
	0		0		0		0
	0  		0		0   	0	
	0  		0		0 		0];	

%----- dMdxpv:
d2BodP2=d2OilFVFdP2(PP,P,T,d60,dg); d2BodPdT=d2OilFVFdPdT(PP,P,T,d60,dg);
d2BodTdP=d2OilFVFdTdP(PP,P,T,d60,dg); d2BodT2=d2OilFVFdT2(PP,P,T,d60,dg);

d2BwdP2=d2WaterFVFdP2(PP,T); d2BwdPdT=d2WaterFVFdPdT(PP,P,T);
d2BwdTdP=d2WaterFVFdTdP(PP,P,T); d2BwdT2=d2WaterFVFdT2(PP,P,T);

d2BgdP2=d2GasFVFdP2(PP,P,T,dg); d2BgdPdT=d2GasFVFdPdT(PP,P,T,dg);
d2BgdTdP=d2GasFVFdTdP(PP,P,T,dg); d2BgdT2=d2GasFVFdT2(PP,P,T,dg);

dMdxpv(:,:,1)=[...	%dMdhL
	2*C*(1/2)/sqrt(hL*(D-hL))*(-2*hL+D)  	    0	0  												    0
	0											0	0												    0
	0  											0	0   											    0
	0  											0	0 												    0];

dMdxpv(:,:,2)=[...	%dMdBSW
	0 											0	0  												    0
	0											0	0												    0
	0  											0	0   											    0
	0  											0	0 												    0];

dMdxpv(:,:,3)=[...	%dMdP
	0 											0	-VWLstd*d2BwdP2-VOLstd*d2BodP2  				    -VWLstd*d2BwdTdP-VOLstd*d2BodTdP
	0											0	0												    0
	0  											0	VGstd*d2BgdP2+VWLstd*d2BwdP2+VOLstd*d2BodP2   	    VGstd*d2BgdTdP+VWLstd*d2BwdTdP+VOLstd*d2BodTdP
	0  											0	0 												    0];



dMdxpv(:,:,4)=[...	%dMdT
	0 											0	-VWLstd*d2BwdPdT-VOLstd*d2BodPdT  					-VWLstd*d2BwdT2-VOLstd*d2BodT2
	0											0	0													0
	0  											0	VGstd*d2BgdPdT+VWLstd*d2BwdPdT+VOLstd*d2BodPdT   	VGstd*d2BgdT2+VWLstd*d2BwdT2+VOLstd*d2BodT2
	0  											0	0 													0];

%----- dNdxcs:
dNdxcs(:,:,1)=[...	%dNdVOLstd
	0 		0		0  		0
	0		0		0		0
	0  		0		0   	0
	0  		0		0 		0];

dNdxcs(:,:,2)=[...	%dNdVWLstd
	0 		0		0  		0
	0		0		0		0
	0  		0		0   	0
	0  		0		0 		0];

dNdxcs(:,:,3)=[...	%dNdVGstd
	0 		0		0  		0
	0		0		0		0
	0  		0		0   	0
	0  		0		0 		0];

dNdxcs(:,:,4)=[...	%dNdE
	0 		0		0  		0
	0		0		0		0
	0  		0		0   	0
	0  		0		0 		0];

%----- dNdxpv:
dNdxpv(:,:,1)=[...	%dNdVhL
	0 						    0				0  		0
	0						    0				0		0
	0  						    0				0   	0
	0  						    0				0 		0];

dNdxpv(:,:,2)=[...	%dNdBSW
	0 							0				0  		0
	-eye(1)						-eye(1)			0		0
	0  							0				0   	0
	0  							0				0 		0];

dNdxpv(:,:,3)=[...	%dNdP
	dBodP 						dBwdP			0  		0
	0							0				0		0
	-dBodP  					-dBwdP			-dBgdP	0
	0  							0				0 		0];

dNdxpv(:,:,4)=[...	%dNdT
	dBodT 						dBwdT			0  		0
	0							0				0		0
	-dBodT  					-dBwdT			-dBgdT	0
	-d60*PP.rho_H2Ostd*cPo  	-rWstd*cPw		0 		0];

%----- dfdx for sys output flag=10
for i=1:4
	dfdx(:,i) = M\(N*dfcsdxpv(:,i) + dNdxpv(:,:,i)*fcs - dMdxpv(:,:,i)*fpv);
	for j=1:4
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
dfcsdu(:,1)=[...	%dfcsdQf_in
	((1-BSW_in))/(24*3600); ...
	(BSW_in)/(24*3600); ...
	((1-BSW_in)*Rs_in)/(24*3600); ...
	((1-BSW_in)*d60*PP.rho_H2Ostd*cPo*(T_in-PP.Tref) + BSW_in*rWstd*cPw*(T_in-PP.Tref))/(24*3600)];

dfcsdu(:,2)=[...	%dfcsdQg_in
	0; ...
	0; ...
	1/(24*3600); ...
	(dg*PP.rho_AIRstd*(eGres_in+DeOG+cPg*(T_in-PP.Tref)))/(24*3600)];
	
dfcsdu(:,3)=[...	%dfcsdBSW_in
	(-Qf_in)/(24*3600); ...
	(Qf_in)/(24*3600); ...
	(-Qf_in*Rs_in)/(24*3600); ...
	(-Qf_in*d60*PP.rho_H2Ostd*cPo*(T_in-PP.Tref) + Qf_in*rWstd*cPw*(T_in-PP.Tref))/(24*3600)];		

dRs_indT_in=dSolGasOilRatiodT(PP,P,T_in,d60,dg);	
deGres_indT_in=dGasSpecResEnthdT(PP,P,T_in,dg);
dfcsdu(:,4)=[...	%dfcsdT_in
	0; ...
	0; ...
	(Qf_in*(1-BSW_in)*dRs_indT_in)/(24*3600); ...
	(Qf_in*(1-BSW_in)*d60*PP.rho_H2Ostd*cPo + Qf_in*BSW_in*rWstd*cPw + Qg_in*dg*PP.rho_AIRstd*(deGres_indT_in+cPg))/(24*3600)];

dfcsdu(:,5)=[...	%dfcsdQf
	(-(1-BSW))/(24*3600); ...
	(-BSW)/(24*3600); ...
	(-(1-BSW)*Rs)/(24*3600); ...
	(-(1-BSW)*d60*PP.rho_H2Ostd*cPo*(T-PP.Tref) - BSW*rWstd*cPw*(T-PP.Tref))/(24*3600)];

dfcsdu(:,6)=[...	%dfcsdQg
	0; ...
	0; ...
	(-1)/(24*3600); ...
	(-dg*PP.rho_AIRstd*(eGres+DeOG+cPg*(T-PP.Tref)))/(24*3600)];

%----- dfdu for sys output flag=11	
dfdu=M\N*dfcsdu;
sys=dfdu;
%--------------------------------------------------------------------
    case 12,	%dhdx
dhdx=eye(4,4);
sys=dhdx;	
%--------------------------------------------------------------------
    case 13,	%dhdu
dhdu=zeros(4,6);
sys=dhdu;
%--------------------------------------------------------------------
    otherwise   % Unexpected flags
DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
end
%--------------------------------------------------------------------    
end