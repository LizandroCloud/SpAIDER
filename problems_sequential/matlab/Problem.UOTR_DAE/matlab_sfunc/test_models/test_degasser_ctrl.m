%Script to run degasser_ctrl Level-1 MATLAB S-Function

%Blackoil object class property package
global PP; PP=blackoil;

%Setting model parameters
d60=141.5/(30+131.5);
ws=(255/PP.rho_H2Ostd)/(1+(255/PP.rho_H2Ostd));
Tbn=(11.8*d60)^3/1.8;
dg=1.05;

par=cell(13,1);
par{1}=8.00;              %C    Vessel length                [m]
par{2}=2.50;              %D    Vessel diameter              [m]
par{3}=d60;               %d60  Oil specific gravity         []
par{4}=Tbn;               %Tbn  Oil normal boiling point     [K]
par{5}=1.05;              %dg   Gas specific gravity         []
par{6}=ws;                %ws   Water salt mass faction      []
par{7}=1800;              %CvL  Liq. valve constant          [gpm]
par{8}=0.91;              %FL   Liq. pressure rec. factor    []
par{9}=0;                 %tauL Valve time constant          [s]
par{10}=120;              %CvG  Gas valve constant           [gpm]
par{11}=0.80;             %Fgm  Spec. heat cap. ratio factor []
par{12}=0.70;             %DPrT Pressure drop ratio factor   []
par{13}=0;                %tauG Gas valve time constant      [s]

%Consistent initial conditions for steady-state case
Qf_in=22000;              %Inlet liquid std. vol. flow rate  [Sm^3/d]
BSW_in=0.15;              %Inlet basic sediments & water     []
Qg_in=5*Qf_in*(1-BSW_in); %Inlet gas std. vol. flow rate     [Sm^3/d]
T_in=80+273.15;           %Inlet temperature                 [K]
P_outL=301.325;           %Liq. valve outlet pressure        [kPa]
P_outG=701.325;           %Gas valve outlet pressure         [kPa]
hL=1.50;                  %Liquid level                      [m]
P=801.325;                %Pressure                          [kPa]
css=[hL; P; Qf_in; Qg_in; BSW_in; T_in; P_outL; P_outG];
[~,x0]=degasser_ctrl([],[],[],0,css,par);

MV_inL=x0(5);
MV_inG=x0(6);

%Getting state-space linear model
u0=[Qf_in; Qg_in; BSW_in; T_in; P_outL; MV_inL; P_outG; MV_inG];
sys=sfunc2ss(@degasser_ctrl,x0,u0,par)
