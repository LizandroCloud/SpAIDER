function [c, ceq]=cnstr(~)
global t x Ndx
% x(5,1:20)
xt=(t(1):(t(end)-t(1))/Ndx:t(end))';

%x=modelo_eq(u);
%Nx=size(x,2);
Nxt=size(xt,1); %Nxt=size(t,1);
ceq=[]; c=[];
%----- Equality constraints ---------------------------------------
%Equality constraint for terminal conditions of states derivatives
ceqdfterm=[...
    ...(t(end)-t(end-1))\(x(end,1:end-1)-x(end-1,1:end-1)),			... from x(1) to x(37)
	(t(end)-t(end-1))\(x(end,end)-x(end-1,end)) - 1000/(3600*24)];	%dVg,flare/dt = 1000 Mm^3/d
ceq=ceqdfterm';
%--------------------------------------------------------------------

%----- Inequality constraints ---------------------------------------
%States upper and lower bounds
fLB = [...
	...0.1;...		%x(1):  SG_001.hWCS   				Water phase level in CS          	[m]
	...0;...		%x(2):  SG_001.BSWLCS				Basic sediments & water in CS 		[]
	...1;...		%x(3):  SG_001.hCS      			Total level in CS                  	[m]
	1;...		%x(4):  SG_001.hCL      			Liquid level in CL                  [m]
	...0;...		%x(5):  SG_001.BSWLCL				Basic sediments & water in CL 		[]
	...1600;...	%x(6):  SG_001.P        			Pressure                            [kPa]
	...270;...		%x(7):  SG_001.T        			Temperature                         [K]
	0.5;...		%x(8):  V_TO_001.hL     			Liquid level                       	[m]
	...0;...		%x(9):  V_TO_001.BSW				Basic sediments & water 			[]
	...700;...		%x(10): V_TO_001.P      			Pressure                            [kPa]
	...340;...		%x(11): V_TO_001.T      			Temperature                         [K]
	...0.1;...		%x(12): TO_001.hW      				Water phase level                   [m]
	...340;...		%x(13): TO_001.T       				Temperature                         [K]
	0.5];...		%x(14): V_TO_002.hL     			Liquid level                       	[m]
	...0;...		%x(15): V_TO_002.BSW				Basic sediments & water 			[]
	...150;...		%x(16): V_TO_002.P      			Pressure                            [kPa]
	...340;...		%x(17): V_TO_002.T      			Temperature                         [K]
	...0.1];...		%x(18): TO_002.hW      				Water phase level                   [m]
	...340;...		%x(19): TO_002.T       				Temperature                         [K]
	...150;...		%x(20): UVRC.scrb(1).P 				Pressure 							[kPa]
	...290;...		%x(21): UVRC.scrb(1).T  			Temperature                         [K]
	...700;...		%x(22): UVRC.scrb(2).P 				Pressure 							[kPa]
	...290;...		%x(23): UVRC.scrb(2).T  			Temperature                         [K]
	...1600;...	%x(35): V_001.P 					Pressure 							[kPa]
	...270];       %x(36): V_001.T  					Temperature                         [K]

fUB = [...
	...0.5; ...	%x(1):  SG_001.hWCS   				Water phase level in CS          	[m]
	...0.5; ...	%x(2):  SG_001.BSWLCS				Basic sediments & water in CS 		[]
	...3; ...		%x(3):  SG_001.hCS      			Total level in CS                  	[m]
	3; ...		%x(4):  SG_001.hCL      			Liquid level in CL                  [m]
	...0.5; ...	%x(5):  SG_001.BSWLCL				Basic sediments & water in CL 		[]
	...2500; ...	%x(6):  SG_001.P        			Pressure                            [kPa]
	...320; ...	%x(7):  SG_001.T        			Temperature                         [K]
	3; ...		%x(8):  V_TO_001.hL     			Liquid level                       	[m]
	...0.5; ...	%x(9):  V_TO_001.BSW				Basic sediments & water 			[]
	...900; ...	%x(10): V_TO_001.P      			Pressure                            [kPa]
	...370; ...	%x(11): V_TO_001.T      			Temperature                         [K]
	...3; ...		%x(12): TO_001.hW      				Water phase level                   [m]
	...370; ...	%x(13): TO_001.T       				Temperature                         [K]
	3]; ...		%x(14): V_TO_002.hL     			Liquid level                       	[m]
	...0.1; ...	%x(15): V_TO_002.BSW				Basic sediments & water 			[]
	...300; ...	%x(16): V_TO_002.P      			Pressure                            [kPa]
	...370; ...	%x(17): V_TO_002.T      			Temperature                         [K]
	...3]; ...		%x(18): TO_002.hW      				Water phase level                   [m]
	...370; ...	%x(19): TO_002.T       				Temperature                         [K]
	...300; ...	%x(20): UVRC.scrb(1).P 				Pressure 							[kPa]
	...330; ...	%x(21): UVRC.scrb(1).T  			Temperature                         [K]
	...900; ...	%x(22): UVRC.scrb(2).P 				Pressure 							[kPa]
	...330; ...	%x(23): UVRC.scrb(2).T  			Temperature                         [K]
	...2500; ...	%x(35): V_001.P 					Pressure 							[kPa]
	...330];       %x(36): V_001.T  					Temperature                         [K]

%idx=[1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;35;36];
idx=[4;8;14];
idxt=zeros(Nxt,1);
for k=1:Nxt
    idxt(k)=max(find(t<=xt(k)));
end
%Lower bounds states' trajectory inequality constraints
cLB=-(x(idxt,idx)-ones(Nxt,1)*fLB');
cLB=cLB(:);

%Upper bounds states' trajectory inequality constraints
cUB=(x(idxt,idx)-ones(Nxt,1)*fUB');
cUB=cUB(:);

%outputting inequality constraints
c=[cLB; cUB];

end
