function [xss,uss]=uotr_ss(x,u,yss,par)
%uotr_ss - Oil treatment unit steady-state condition computation
%
%This function updates the values in x and u after solving the 
%steady-state condition for the oil treatment unit upstream topside 
%facility model, described by the uotr.m S-Function file.
%
%x and u are, respectively, the state and input variables passed from
%the user, which are initial guess for the steady-state computation.
%yss is the steady-state condition for the uotr.m outputs identified
%as controlled variables after building a multi-loop control scheme.
%par are the uotr.m model parameters.
%
%xss and uss are, respectively, the updated state and input variables
%at steady-state condition. The variables in xss updated are the ones
%which remain free after building the multi-loop control scheme. The 
%variables in uss updated are manipulated variables of the multi-loop
%control scheme.
%
%Fixed states: 	Controlled variables. Values not updated.
%Fixed inputs:	Exogenous input variables. Values not updated.
%Fixed outpus:  Controlled variables. Steady-state condition as yss.
%Free states:  	Free variables. Values updated.
%Free inputs: 	Manipulated variables. Values updated.
%
%----- States -------------------------------------------------------
%Fixed:
%x(1):  SG_001.hWCS   				Water phase level in CS          		[m]
%x(4):  SG_001.hCL      			Liquid level in CL                  	[m]
%x(8):  V_TO_001.hL     			Liquid level                       		[m]
%x(11): V_TO_001.T      			Temperature                         	[K]
%x(12): TO_001.hW      				Water phase level                   	[m]
%x(14): V_TO_002.hL     			Liquid level                       		[m]
%x(16): V_TO_002.P      			Pressure                            	[kPa]
%x(18): TO_002.hW      				Water phase level                   	[m]
%x(21): UVRC.scrb(1).T  			Temperature                         	[K]
%x(23): UVRC.scrb(2).T  			Temperature                         	[K]
%x(26): V_001.P 					Pressure 								[kPa]
%x(27): V_001.T  					Temperature                         	[K]
%x(28): Vg,flare       				Totalized flaring gas std. vol.     	[SMm3]

%Free:
%x(2):  SG_001.BSWLCS				Basic sediments & water in CS 			[]
%x(3):  SG_001.hCS      			Total level in CS                  		[m]
%x(5):  SG_001.BSWLCL				Basic sediments & water in CL 			[]
%x(6):  SG_001.P        			Pressure                            	[kPa]
%x(7):  SG_001.T        			Temperature                         	[K]
%x(9):  V_TO_001.BSW				Basic sediments & water 				[]
%x(10): V_TO_001.P      			Pressure                            	[kPa]
%x(13): TO_001.T       				Temperature                         	[K]
%x(15): V_TO_002.BSW				Basic sediments & water 				[]
%x(17): V_TO_002.T      			Temperature                         	[K]
%x(19): TO_002.T       				Temperature                         	[K]
%x(20): UVRC.scrb(1).P 				Pressure 								[kPa]
%x(22): UVRC.scrb(2).P 				Pressure 								[kPa]
%x(24): UVRC.asc(1).MVinput     	MV tracking       						[]
%x(25): UVRC.asc(2).MVinput     	MV tracking       						[]
%--------------------------------------------------------------------
%----- Inputs -------------------------------------------------------
%Fixed:
%u(1):  feed.Qf_in     		Liq. std. vol. flow rate    		[Sm^3/d]
%u(2):  feed.Rgo     		Gas-oil ratio			       		[Sm^3/Sm^3]
%u(3):  feed.BSW_in    		Basic sediments & water       		[]
%u(4):  feed.T_in      		Temperature                   		[K]
%u(6):  LV_002.P_out    	Outlet pressure            			[kPa]
%u(8):  PV_001.P_out    	Outlet pressure     				[kPa]
%u(10): stk.P    			Stock tank pressure            		[kPa]
%u(11): P_001.DP_C    		Cold side pressure drop				[kPa]
%u(12): P_001.DP_H    		Hot side pressure drop				[kPa]
%u(13): P_002.DP    		Pressure drop						[kPa]
%u(17): dslt.Qw    	 		Desalt water liq. std. vol. flow 	[Sm^3/d]
%u(18): dslt.T    	 		Desalt water temperature			[K]
%u(22):	UVRC.cool(1).DP   	Pressure drop        				[kPa]
%u(24):	UVRC.cool(2).DP   	Pressure drop        				[kPa]
%u(27):	P_003.DP     		Pressure drop 						[kPa]

%Free:
%u(5):  LV_001.c    		Valve stem position   				[]
%u(7):  LV_002.c    		Valve stem position 				[]
%u(9):  PV_001.c    		Stem position						[]
%u(14): P_002.Q_in      	Heat applied to the system          [kW]
%u(15): LV_003.c    		Valve stem position   				[]
%u(16): LV_004.c    		Valve stem position 				[]
%u(19): LV_005.c    		Valve stem position   				[]
%u(20): LV_006.c    		Valve stem position 				[]
%u(21):	UVRC.cool(1).Q_out 	Heat removed from the system        [kW]
%u(23):	UVRC.cool(2).Q_out  Heat removed from the system        [kW]
%u(25):	UVRC.MV_in     		Input signal from controller        []
%u(26):	P_003.Q_out     	Heat removed from the system        [kW]
%--------------------------------------------------------------------
%------ Outputs for steady-state condition --------------------------
%Fixed:
%There are no fixed output variables
%--------------------------------------------------------------------

%2023. Version 0. Coded by Diego Trica (diegotrica@gmail.com)

%----- Fixed and free states:
idx = [...
	true; ...	 SG_001.hWCS
	false; ...	 SG_001.BSWLCS
	false; ...	 SG_001.hCS
	true; ...	 SG_001.hCL
	false; ...	 SG_001.BSWLCL
	false; ...	 SG_001.P
	false; ...	 SG_001.T
	true; ...	 V_TO_001.hL
	false; ...	 V_TO_001.BSW
	false; ...	 V_TO_001.P
	true; ...	 V_TO_001.T
	true; ...	 TO_001.hW
	false; ... 	 TO_001.T
	true; ...	 V_TO_002.hL
	false; ...	 V_TO_002.BSW
	true; ...	 V_TO_002.P
	false; ...	 V_TO_002.T
	true; ...	 TO_002.hW
	false; ...	 TO_002.T
	false; ...	 UVRC.scrb(1).P
	true; ...	 UVRC.scrb(1).T
	false; ...	 UVRC.scrb(2).P
	true; ...	 UVRC.scrb(2).T
	false; ... 	 UVRC.asc(1).MVinput
	false; ... 	 UVRC.asc(2).MVinput
	true; ...	 V_001.P
	true ...	 V_001.T
%	true ...	 Vg,flare
	];
fixed_x	= x(idx==true);
free_x 	= x(idx==false);

%----- Fixed and free inputs:
idu = [...
	true; ...	feed.Qf_in  
	true; ...	feed.Rgo 
	true; ...	feed.BSW_in
	true; ...	feed.T_in  
	false; ...	LV_001.c
	true; ...	LV_002.P_out
	false; ...	LV_002.c
	true; ...	PV_001.P_out
	false; ...	PV_001.c
	true; ...	stk.P
	true; ...	P_001.DP_C
	true; ...	P_001.DP_H
	true; ...	P_002.DP
	false; ...	P_002.Q_in
	false; ...	LV_003.c
	false; ...	LV_004.c
	true; ...	dslt.Qw
	true; ...	dslt.T
	false; ...	LV_005.c
	false; ...	LV_006.c
	false; ...	UVRC.cool(1).Q_out
	true; ...	UVRC.cool(1).DP
	false; ...	UVRC.cool(2).Q_out
	true; ...	UVRC.cool(2).DP
	false; ...	UVRC.MV_in
	false; ...	P_003.Q_out
	true ...	P_003.DP
	];
fixed_u = u(idu==true);
free_u 	= u(idu==false);

%----- Fixed outputs:
idy 	= zeros(37,1);

%Steady-state computation:
%Find feasible region
opt = optimoptions('fsolve','Algorithm','trust-region-reflective','Display','iter','MaxIter',100,'TolFun',1e-2);
z = fsolve(...
	@(z)fh(z,fixed_x,idx,fixed_u,idu,yss,idy,par), ...
	[free_x;free_u],opt);
if imag(z)>1e-4;
	error('solution is imag')
end
free_x = real(z(1:size(find(idx==false),1)));
free_u = real(z(size(find(idx==false),1)+1:end));

%Refine to find solution
opt = optimoptions('fsolve','Algorithm','trust-region-dogleg','Display','iter','MaxFunEvals',25000,'MaxIter',500);
z = fsolve(...
	@(z)fh(z,fixed_x,idx,fixed_u,idu,yss,idy,par), ...
	[free_x;free_u],opt);
if imag(z)>1e-4;
	error('solution is imag')
end
free_x = real(z(1:size(find(idx==false),1)));
free_u = real(z(size(find(idx==false),1)+1:end));

%Allocating solution
xss = x; uss = u;
xss(idx==false) = free_x;
uss(idu==false) = free_u;
end

function res = fh(z,fixed_x,idx,fixed_u,idu,yss,idy,par)
    free_x          = z(1:size(find(idx==false),1));
    free_u          = z(size(find(idx==false),1)+1:end);

    x(idx==true) 	= fixed_x;
    x(idx==false) 	= free_x;
	x = x(:);
	
    u(idu==true) 	= fixed_u;
	u(idu==false) 	= free_u;
	u = u(:);
	
	res = uotr([],x,u,1,[],par);
	res = res(1:end-1);
	if ~isempty(yss)
		y = uotr([],x,u,3,[],par);
		res = [res;(yss-y(idy==true))];
	end
end