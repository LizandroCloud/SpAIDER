function [ceq c dceq dc]=cr2(f,t_ult,t_out,u,y,stages)
%global t x Ndx
%global x

%xt=(t(1):(t(end)-t(1))/Ndx:t(end))';
Nut=size(ut,1); Nxt=size(xt,1);
[x,G] = modelo_eq(t,u,ut,ufinal,modelo_x0,modelo_u0,modelo_dFdx,modelo_dFdu);
c=[]; ceq=[]; gradc=[]; gradceq=[]; 
%----- Equality constraints ---------------------------------------
%Equality constraint for terminal conditions of states derivatives

ceqdfterm=[...
    ...(t(end)-t(end-1))\(x(end,1:end-1)-x(end-1,1:end-1)),			... from x(1) to x(37)
	-(t(end)-t(end-1))\(x(end,end)-x(end-1,end)) - 1000/(3600*24)];	%dVg,flare/dt = 1000 Mm^3/d
ceq=ceqdfterm';

%Gradient
g=zeros(Nut-1,1); gh=zeros(Nut-1,1);
for m=1:Nut-1
	idut=max(find(t<=ut(m))); iduth=max(find(t(1:end-1)<=ut(m)));
	idutp=max(find(t<=ut(m+1))); idutph=max(find(t(1:end-1)<=ut(m+1)));
	for n=idut:idutp-1
		g(m) = g(m) + G(end,1,n);
	end
	for n=iduth:idutph-1
		gh(m) = gh(m) + G(end,1,n);
	end	
end
gradceq=-(t(end)-t(end-1))\(g-gh)*1e3; %Retornando ordem de grandeza

%--------------------------------------------------------------------

%----- Inequality constraints ---------------------------------------
%{
%States upper and lower bounds
fLB = [...
	0.1;...		%x(1):  SG_001.hWCS   				Water phase level in CS          	[m]
	0;...		%x(2):  SG_001.BSWLCS				Basic sediments & water in CS 		[]
	1;...		%x(3):  SG_001.hCS      			Total level in CS                  	[m]
	1;...		%x(4):  SG_001.hCL      			Liquid level in CL                  [m]
	0;...		%x(5):  SG_001.BSWLCL				Basic sediments & water in CL 		[]
	1600;...	%x(6):  SG_001.P        			Pressure                            [kPa]
	270;...		%x(7):  SG_001.T        			Temperature                         [K]
	0.5;...		%x(8):  V_TO_001.hL     			Liquid level                       	[m]
	0;...		%x(9):  V_TO_001.BSW				Basic sediments & water 			[]
	700;...		%x(10): V_TO_001.P      			Pressure                            [kPa]
	340;...		%x(11): V_TO_001.T      			Temperature                         [K]
	0.1;...		%x(12): TO_001.hW      				Water phase level                   [m]
	340;...		%x(13): TO_001.T       				Temperature                         [K]
	0.5;...		%x(14): V_TO_002.hL     			Liquid level                       	[m]
	0;...		%x(15): V_TO_002.BSW				Basic sediments & water 			[]
	150;...		%x(16): V_TO_002.P      			Pressure                            [kPa]
	340;...		%x(17): V_TO_002.T      			Temperature                         [K]
	0.1;...		%x(18): TO_002.hW      				Water phase level                   [m]
	340;...		%x(19): TO_002.T       				Temperature                         [K]
	150;...		%x(20): UVRC.scrb(1).P 				Pressure 							[kPa]
	290;...		%x(21): UVRC.scrb(1).T  			Temperature                         [K]
	700;...		%x(22): UVRC.scrb(2).P 				Pressure 							[kPa]
	290;...		%x(23): UVRC.scrb(2).T  			Temperature                         [K]
	1600;...	%x(27): V_001.P 					Pressure 							[kPa]
	270];       %x(28): V_001.T  					Temperature                         [K]

fUB = [...
	0.5; ...	%x(1):  SG_001.hWCS   				Water phase level in CS          	[m]
	0.5; ...	%x(2):  SG_001.BSWLCS				Basic sediments & water in CS 		[]
	3; ...		%x(3):  SG_001.hCS      			Total level in CS                  	[m]
	3; ...		%x(4):  SG_001.hCL      			Liquid level in CL                  [m]
	0.5; ...	%x(5):  SG_001.BSWLCL				Basic sediments & water in CL 		[]
	2500; ...	%x(6):  SG_001.P        			Pressure                            [kPa]
	320; ...	%x(7):  SG_001.T        			Temperature                         [K]
	3; ...		%x(8):  V_TO_001.hL     			Liquid level                       	[m]
	0.5; ...	%x(9):  V_TO_001.BSW				Basic sediments & water 			[]
	900; ...	%x(10): V_TO_001.P      			Pressure                            [kPa]
	370; ...	%x(11): V_TO_001.T      			Temperature                         [K]
	3; ...		%x(12): TO_001.hW      				Water phase level                   [m]
	370; ...	%x(13): TO_001.T       				Temperature                         [K]
	3; ...		%x(14): V_TO_002.hL     			Liquid level                       	[m]
	0.1; ...	%x(15): V_TO_002.BSW				Basic sediments & water 			[]
	320; ...	%x(16): V_TO_002.P      			Pressure                            [kPa]
	370; ...	%x(17): V_TO_002.T      			Temperature                         [K]
	3; ...		%x(18): TO_002.hW      				Water phase level                   [m]
	370; ...	%x(19): TO_002.T       				Temperature                         [K]
	320; ...	%x(20): UVRC.scrb(1).P 				Pressure 							[kPa]
	330; ...	%x(21): UVRC.scrb(1).T  			Temperature                         [K]
	900; ...	%x(22): UVRC.scrb(2).P 				Pressure 							[kPa]
	330; ...	%x(23): UVRC.scrb(2).T  			Temperature                         [K]
	2500; ...	%x(27): V_001.P 					Pressure 							[kPa]
	330];       %x(28): V_001.T  					Temperature                         [K]
%}
%States upper and lower bounds
fLB = [...
	0.1;...		%x(1):  SG_001.hWCS   				Water phase level in CS          	[m]
	0;...		%x(2):  SG_001.BSWLCS				Basic sediments & water in CS 		[]
	2.05;...	%x(3):  SG_001.hCS      			Total level in CS                  	[m]
	1.2;...		%x(4):  SG_001.hCL      			Liquid level in CL                  [m]
	0;...		%x(5):  SG_001.BSWLCL				Basic sediments & water in CL 		[]
	1950;...	%x(6):  SG_001.P        			Pressure                            [kPa]
	300;...		%x(7):  SG_001.T        			Temperature                         [K]
	0.90;...	%x(8):  V_TO_001.hL     			Liquid level                       	[m]
	0;...		%x(9):  V_TO_001.BSW				Basic sediments & water 			[]
	750;...		%x(10): V_TO_001.P      			Pressure                            [kPa]
	350;...		%x(11): V_TO_001.T      			Temperature                         [K]
	1.8;...		%x(12): TO_001.hW      				Water phase level                   [m]
	350;...		%x(13): TO_001.T       				Temperature                         [K]
	0.90;...	%x(14): V_TO_002.hL     			Liquid level                       	[m]
	0;...		%x(15): V_TO_002.BSW				Basic sediments & water 			[]
	250;...		%x(16): V_TO_002.P      			Pressure                            [kPa]
	350;...		%x(17): V_TO_002.T      			Temperature                         [K]
	1.8;...		%x(18): TO_002.hW      				Water phase level                   [m]
	350;...		%x(19): TO_002.T       				Temperature                         [K]
	250;...		%x(20): UVRC.scrb(1).P 				Pressure 							[kPa]
	290;...		%x(21): UVRC.scrb(1).T  			Temperature                         [K]
	750;...		%x(22): UVRC.scrb(2).P 				Pressure 							[kPa]
	290;...		%x(23): UVRC.scrb(2).T  			Temperature                         [K]
	1930;...	%x(27): V_001.P 					Pressure 							[kPa]
	300];       %x(28): V_001.T  					Temperature                         [K]

fUB = [...
	0.3; ...	%x(1):  SG_001.hWCS   				Water phase level in CS          	[m]
	0.5; ...	%x(2):  SG_001.BSWLCS				Basic sediments & water in CS 		[]
	2.25; ...	%x(3):  SG_001.hCS      			Total level in CS                  	[m]
	1.8; ...	%x(4):  SG_001.hCL      			Liquid level in CL                  [m]
	0.5; ...	%x(5):  SG_001.BSWLCL				Basic sediments & water in CL 		[]
	2150; ...	%x(6):  SG_001.P        			Pressure                            [kPa]
	315; ...	%x(7):  SG_001.T        			Temperature                         [K]
	1.1; ...	%x(8):  V_TO_001.hL     			Liquid level                       	[m]
	0.5; ...	%x(9):  V_TO_001.BSW				Basic sediments & water 			[]
	850; ...	%x(10): V_TO_001.P      			Pressure                            [kPa]
	370; ...	%x(11): V_TO_001.T      			Temperature                         [K]
	2.2; ...	%x(12): TO_001.hW      				Water phase level                   [m]
	370; ...	%x(13): TO_001.T       				Temperature                         [K]
	1.1; ...	%x(14): V_TO_002.hL     			Liquid level                       	[m]
	0.1; ...	%x(15): V_TO_002.BSW				Basic sediments & water 			[]
	320; ...	%x(16): V_TO_002.P      			Pressure                            [kPa]
	370; ...	%x(17): V_TO_002.T      			Temperature                         [K]
	2.2; ...	%x(18): TO_002.hW      				Water phase level                   [m]
	370; ...	%x(19): TO_002.T       				Temperature                         [K]
	320; ...	%x(20): UVRC.scrb(1).P 				Pressure 							[kPa]
	325; ...	%x(21): UVRC.scrb(1).T  			Temperature                         [K]
	850; ...	%x(22): UVRC.scrb(2).P 				Pressure 							[kPa]
	325; ...	%x(23): UVRC.scrb(2).T  			Temperature                         [K]
	2150; ...	%x(27): V_001.P 					Pressure 							[kPa]
	315];       %x(28): V_001.T  					Temperature                         [K]

idx=[1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;27;28];
%idx=[4;6;8;10;12;14;16;18;20;22;27];
idxt=zeros(Nxt,1);
for k=1:Nxt
    idxt(k)=max(find(t<=xt(k)));
end

%Lower bounds states' trajectory inequality constraints
%cLB=-(x(idxt,idx)-ones(Nxt,1)*fLB');
%cLB=cLB(:);

%cLB=[]; gradcLB=[];
%for j=1:length(idx)
%	for k=1:Nxt
%		cLB=[cLB;-(x(idxt(k),idx(j))-fLB(j))];
%        g=G(idx(j),1,1:Nut-1); g=g(:);
%        gradcLB=[gradcLB,-g];
%    end
%end

cLB=[]; gradcLB=[];
for j=1:length(idx)
	for k=1:Nxt
		cLB=[cLB;-(x(idxt(k),idx(j))-fLB(j))];
        g=zeros(Nut-1,1);
		for m=1:Nut-1
			idut=max(find(t(1:idxt(k))<=ut(m)));
			idutp=max(find(t(1:idxt(k))<=ut(m+1)));
			for n=idut:idutp-1
				g(m) = g(m) + G(idx(j),1,n);
			end
		end
        gradcLB=[gradcLB,-g];
    end
end	

%Upper bounds states' trajectory inequality constraints
cUB=(x(idxt,idx)-ones(Nxt,1)*fUB');
cUB=cUB(:);
gradcUB=-gradcLB;
%cUB=[]; gradcUB=[];
%for j=1:length(idx)
%	for k=1:Nxt
%		cUB=[cUB;(x(idxt(k),idx(j))-fUB(j))];
%	end
%	g=G(idx(j),1,1:Nut-1); g=g(:)';
%	gradcUB=[gradcUB;g];	
%end

%outputting inequality constraints and gradients
c=[cLB; cUB];
gradc=[gradcLB, gradcUB]*1e3; %Retornando ordem de grandeza

end