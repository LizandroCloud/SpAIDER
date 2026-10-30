function [ dx ] = modelo_eq( t_out, x, ws, ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,tal,est,yi,utarg,tstarg,nstarg)
[u] = ver_reg_est1(ws,yi,t_out);
[uref] = ver_reg_est1(ws,yi,t_out);

X=x;
T = u(1);
R=8.314;
%Reactions' pre-exponential factors
k0(1)=1.8E-4; %[kmol.bar0.5/kgcat.h]
k0(2)=7.6;    %[kmol/bar.kgcat.h]
k0(3)=2.2E-5; %[kmol.bar0.5/kgcat.h]
%Reactions' activation energies [J/mol]
Ea(1)=240.1E3;
Ea(2)=67.1E3;
Ea(3)=243.9E3;
%Reactions' effectiveness factors [-]
eta(1)=0.01;
eta(2)=0.3;
eta(3)=0.01;
%Adsorption constants' pre-exponential factors
Ka0(1)=1.8E-1; %at 823 K [1/bar]
Ka0(2)=0.4;    %at 823 K [-]
Ka0(3)=40.9;   %at 648 K [1/bar]
Ka0(5)=2.9E-2; %at 648 K [1/bar]
%Adsorption constants' reference temperatures
Tref=648*ones(5,1);
Tref(1)=823;
Tref(2)=823;
%Adsorption enthalpies [J/mol]
dH(1)=-38.3E3;
dH(2)=88.7E3;
dH(3)=-70.7E3;
dH(5)=-82.9E3;
%Components' Gibbs free energy constants
Gc=1E3*[-75.3   7.6E-2 1.9E-5   %[G(j,1)]=J/mol
        -241.7  4.2E-2 7.4E-6   %[G(j,2)]=J/mol.K
        -109.9 -9.2E-2 1.5E-6   %[G(j,3)]=J/mol.K2
        -393.4  3.8E-3 1.3E-6];
%Stoichiometric coefficients [-]
mu=[-1  0 -1
    -1 -1 -2
     1 -1  0
     0  1  1
     3  1  4];

%Reactor's physical parameters
%-----------------------------------------------------------------------------
L=7;                   %Reactor length [m]
r0=2*2.54E-2;          %Permeation zone's outer radius [m]
delta=5*1E-6;          %Membrane thickness [m]
rit=r0+delta;          %Permeation zone's total radius [m]
Ri=2.5*2.54E-2;        %Reactor radius [m]
Ac=3.1416*(Ri^2-rit^2);   %Reaction zone's cross section area [m2]
epsib=0.5;             %Bed void fraction [-]
rhocat=2100;           %Catalyst density [kg/m3]
rhob=rhocat*(1-epsib); %Bulk density [kg/m3]
%-----------------------------------------------------------------------------
%Operating conditions 
%-----------------------------------------------------------------------------
Pp=1;                  %Permeate total pressure [bar]
Fa0=1;                  % CH4 initial feed, [kmol/h]
FCH40=1;               %Inlet CH? flow [kmol/h]
teta=zeros(5,1) ;       %Inlet ratios [-]
teta(2)=10;             %FH?O0/FCH?0 [-] % VARIÁVEL DE DECSÃO 

teta(5)=0.01;          %FH?0/FCH?0 [-]
sweep=10;              %Fsweep/FCH?0 [-]
Pi=3.14159;
PhiB=0.5; % Bed porosity
%-----------------------------------------------------------------------------
%Integration parameters
%-----------------------------------------------------------------------------
%Boundary conditions
ksi0=0; X0MR=zeros(3,1); X0PBR=zeros(2,1);

%-----------------------------------------------------------------------------

Pt=20; % If Pt is the parameter/Matrix
 
% [k,Ka,Ke]=kin(T)
% other parameters
Ac=Pi*(Ri^2-rit^2); % area, [m2]
Qpd=6.3227*(1e-3)*exp(-15630/(R*T)); % H2 permeability in Pd, [m3.um.h.atm0.5/m2]
rhoB=rhocat*(1-PhiB); % Bed density [kg/m3]
C1=rhoB*Ac*L/Fa0; % Constant for diff. eq. 1 and 2, PBR and MR
C2=Qpd*2*Pi*rit*L/(Fa0*delta*22.4); % Constant for diff. eq. 3, MR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for i=1:3
        k(i)=k0(i)*exp((Ea(i)/R)*(1/648-1/T));
    end
    for j=1:5
        Ka(j)=Ka0(j)*exp((dH(j)/R)*(1/Tref(j)-1/T));
    end
    for j=1:4
        G0(j)=Gc(j,1)+Gc(j,2)*T+Gc(j,3)*T^2;
    end
    dGr=zeros(3,1);
    for i=1:3
        for j=1:4
            dGr(i)=dGr(i)+mu(j,i)*G0(j);
        end
        Ke(i)=exp(-dGr(i)/(R*T));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




    S=Pt/(1+teta(2)+teta(3)+teta(4)+teta(5)+2*X(1)-X(3));
    P(1)=(1-X(1))*S;
    P(2)=(teta(2)-X(1)-X(2))*S;
    P(3)=(teta(3)+X(1)-X(2))*S;
    P(4)=(teta(4)+X(2))*S;
    P(5)=(teta(5)+3*X(1)+X(2)-X(3))*S;
    DEN2=(1+Ka(3)*P(3)+Ka(5)*P(5)+Ka(1)*P(1)+Ka(2)*P(2)/P(5))^2;
    r(1)=eta(1)*(k(1)/(P(5)^2.5))*(P(1)*P(2)-(P(5)^3*P(3))/Ke(1))/DEN2;
    r(2)=eta(2)*(k(2)/P(5))*(P(3)*P(2)-(P(5)*P(4))/Ke(2))/DEN2;
    r(3)=eta(3)*(k(3)/(P(5)^3.5))*(P(1)*P(2)^2-(P(5)^4*P(4))/Ke(3))/DEN2;
    % Pressure of hydrogen in the permeate side
    Pperm=Pp*X(3)/(X(3)+sweep);
    % Differential system
    dx(1)=C1*(r(1)+r(3));          %dXCH?/dksi
    dx(2)=C1*(r(2)+r(3));          %dXCO?/dksi
    dx(3)=C2*(P(5)^0.5-Pperm^0.5); %dYH?/dksi
dx=transp(dx);
end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
