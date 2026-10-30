%Script to run mixer_W3ph Level-1 MATLAB S-Function
clear

%Physical properties
global PP
PP=blackoil;

%Setting model parameters
par=cell(6,1);
par{1}=141.5/(30+131.5);          %d60, Oil specific gravity []
par{2}=(1/1.8)*(11.8*par{1})^3;   %Tbn, Oil normal boiling point [K]
par{3}=255;     				  %Cws, Water salt concentration [g/L]
par{4}=0.7;                       %dg, Gas specific gravity []
par{5}=1;                         %nF, Number of F inlet streams []
par{6}=2;                      	  %nW, Number of W inlet streams []

%Model inputs linearization point
%Inlet F(1):
F(1).Qf		= 22000;			%Std. vol. flow rate 		[Sm^3/d]
F(1).BSW	= 0.005;			%Basic sediments & water 	[]
F(1).Qg		= 200*F(1).Qf*(1-F(1).BSW);
F(1).T     	= 290;				%Temperauture 				[K]

%Inlet W(1)
W(1).Qw		= 5400;				%Std. vol. flow rate 		[Sm^3/d]
W(1).T     	= 360; 				%Temperauture 				[K]

%Inlet W(2)
W(2).Qw		= 1200;				%Std. vol. flow rate 		[Sm^3/d]
W(2).T     	= 360;				%Temperauture 				[K]

P = 2101.325;					%Mixer node pressure 	[kPa]

%Getting state-space linear model
u0 = zeros(par{5}*4+par{6}*2+1,1);

for k=1:par{5}
u0(4*(k-1)+1) = F(k).Qf;
u0(4*(k-1)+2) = F(k).Qg;
u0(4*(k-1)+3) = F(k).BSW;
u0(4*(k-1)+4) = F(k).T;
end

for k=1:par{6}
u0(4*par{5}+2*(k-1)+1) = W(k).Qw;
u0(4*par{5}+2*(k-1)+2) = W(k).T;
end

u0(end) = P;

sys=sfunc2ss(@mixer_W3ph,[],u0,par)