function [xss,uss]=uotr_ctrl_ss(x,u,yss,par)
%uotr_ctrl_ss - %uotr_ctrl.m model steady-state condition computation
%
%This function updates the values in x and u after solving the 
%steady-state condition for the uotr_ctrl.m model S-Function file.
%
%x and u are, respectively, the state and input variables passed from
%the user, which are initial guess for the steady-state computation.
%yss is the steady-state condition for the uotr_ctrl.m outputs iden-
%tified as controlled variables after building a multi-loop control 
%scheme. par are the uotr.m model parameters.
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
%x(1):  SG_001.hWCS   				Water phase level in CS          	[m]
%x(4):  SG_001.hCL      			Liquid level in CL                  [m]
%x(8):  V_TO_001.hL     			Liquid level                       	[m]
%x(11): V_TO_001.T      			Temperature                         [K]
%x(12): TO_001.hW      				Water phase level                   [m]
%x(14): V_TO_002.hL     			Liquid level                       	[m]
%x(16): V_TO_002.P      			Pressure                            [kPa]
%x(18): TO_002.hW      				Water phase level                   [m]
%x(21): UVRC.scrb(1).T  			Temperature                         [K]
%x(23): UVRC.scrb(2).T  			Temperature                         [K]
%x(27): UVRC.TIC(1).derivFilt   	Derivative action with filter      	[]
%x(30): UVRC.TIC(2).derivFilt   	Derivative action with filter      	[]
%x(33): UVRC.MCC.derivFilt   		Derivative action with filter      	[]
%x(35): V_001.P 					Pressure 							[kPa]
%x(36): V_001.T  					Temperature                         [K]
%x(38): LIC_001.derivFilt    		Derivative action with filter       []
%x(41): LIC_002.derivFilt    		Derivative action with filter       []
%x(44): LIC_003.derivFilt    		Derivative action with filter       []
%x(47): LIC_004.derivFilt    		Derivative action with filter       []
%x(50): LIC_005.derivFilt    		Derivative action with filter       []
%x(53): LIC_006.derivFilt    		Derivative action with filter       []
%x(56): PIC_001.derivFilt    		Derivative action with filter       []
%x(59): TIC_001.derivFilt   		Derivative action with filter      	[]
%x(62): TIC_002.derivFilt   		Derivative action with filter      	[]
%x(64): Vg,flare       				Totalized flaring gas std. vol.     [SMm3]

%Free:
%x(2):  SG_001.BSWLCS				Basic sediments & water in CS 		[]
%x(3):  SG_001.hCS      			Total level in CS                  	[m]
%x(5):  SG_001.BSWLCL				Basic sediments & water in CL 		[]
%x(6):  SG_001.P        			Pressure                            [kPa]
%x(7):  SG_001.T        			Temperature                         [K]
%x(9):  V_TO_001.BSW				Basic sediments & water 			[]
%x(10): V_TO_001.P      			Pressure                            [kPa]
%x(13): TO_001.T       				Temperature                         [K]
%x(15): V_TO_002.BSW				Basic sediments & water 			[]
%x(17): V_TO_002.T      			Temperature                         [K]
%x(19): TO_002.T       				Temperature                         [K]
%x(20): UVRC.scrb(1).P 				Pressure 							[kPa]
%x(22): UVRC.scrb(2).P 				Pressure 							[kPa]
%x(24): UVRC.asc(1).MVinput     	MV tracking       					[]
%x(25): UVRC.asc(2).MVinput     	MV tracking       					[]
%x(26): UVRC.TIC(1).AWBTintterm 	Integral action with AWBT  			[]
%x(28): UVRC.TIC(1).MVinput     	MV tracking       					[]
%x(29): UVRC.TIC(2).AWBTintterm 	Integral action with AWBT  			[]
%x(31): UVRC.TIC(2).MVinput     	MV tracking       					[]
%x(32): UVRC.MCC.AWBTintterm 		Integral action with AWBT  			[]
%x(34): UVRC.MCC.MVinput     		MV tracking       					[]
%x(37): LIC_001.AWBTintterm  		Integral action with AWBT  			[]
%x(39): LIC_001.MVinput     		MV tracking       					[]
%x(40): LIC_002.AWBTintterm  		Integral action with AWBT  			[]
%x(42): LIC_002.MVinput     		MV tracking       					[]
%x(43): LIC_003.AWBTintterm  		Integral action with AWBT  			[]
%x(45): LIC_003.MVinput     		MV tracking       					[]
%x(46): LIC_004.AWBTintterm  		Integral action with AWBT  			[]
%x(48): LIC_004.MVinput     		MV tracking       					[]
%x(49): LIC_005.AWBTintterm  		Integral action with AWBT  			[]
%x(51): LIC_005.MVinput     		MV tracking       					[]
%x(52): LIC_006.AWBTintterm  		Integral action with AWBT  			[]
%x(54): LIC_006.MVinput     		MV tracking       					[]
%x(55): PIC_001.AWBTintterm  		Integral action with AWBT  			[]
%x(57): PIC_001.MVinput     		MV tracking       					[]
%x(58): TIC_001.AWBTintterm 		Integral action with AWBT  			[]
%x(60): TIC_001.MVinput     		MV tracking       					[]
%x(61): TIC_002.AWBTintterm 		Integral action with AWBT  			[]
%x(63): TIC_002.MVinput     		MV tracking       					[]
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

%--------------------------------------------------------------------
%------ Outputs for steady-state condition --------------------------
%Fixed:
%There are no fixed output variables
%--------------------------------------------------------------------

%2024-06-10. Version 0. Coded by Diego Trica (diegotrica@gmail.com)

%----- Fixed and free states:
idx = [...
	true; ...	SG_001.hWCS
	false; ...	SG_001.BSWLCS
	false; ...	SG_001.hCS
	true; ...	SG_001.hCL
	false; ...	SG_001.BSWLCL
	false; ...	SG_001.P
	false; ...	SG_001.T
	true; ...	V_TO_001.hL
	false; ...	V_TO_001.BSW
	false; ...	V_TO_001.P
	true; ...	V_TO_001.T
	true; ...	TO_001.hW
	false; ... 	TO_001.T
	true; ...	V_TO_002.hL
	false; ...	V_TO_002.BSW
	true; ...	V_TO_002.P
	false; ...	V_TO_002.T
	true; ...	TO_002.hW
	false; ...	TO_002.T
	false; ...	UVRC.scrb(1).P
	true; ...	UVRC.scrb(1).T
	false; ...	UVRC.scrb(2).P
	true; ...	UVRC.scrb(2).T
	false; ... 	UVRC.asc(1).MVinput
	false; ... 	UVRC.asc(2).MVinput
	false;... 	UVRC.TIC(1).AWBTintterm 
	true;...	UVRC.TIC(1).derivFilt
	false;...	UVRC.TIC(1).MVinput
	false;... 	UVRC.TIC(2).AWBTintterm
	true;...	UVRC.TIC(2).derivFilt
	false;...	UVRC.TIC(2).MVinput
	false;... 	UVRC.MCC.AWBTintterm
	true;...	UVRC.MCC.derivFilt
	false;...	UVRC.MCC.MVinput
	true; ...	V_001.P
	true; ...	V_001.T
	false;... 	LIC_001.AWBTintterm 
	true;...	LIC_001.derivFilt   
	false;...	LIC_001.MVinput     
	false;... 	LIC_002.AWBTintterm 
	true;...	LIC_002.derivFilt   
	false;...	LIC_002.MVinput     
	false;... 	LIC_003.AWBTintterm 
	true;...	LIC_003.derivFilt   
	false;...	LIC_003.MVinput     
	false;... 	LIC_004.AWBTintterm 
	true;...	LIC_004.derivFilt   
	false;...	LIC_004.MVinput     
	false;... 	LIC_005.AWBTintterm 
	true;...	LIC_005.derivFilt   
	false;...	LIC_005.MVinput     
	false;... 	LIC_006.AWBTintterm 
	true;...	LIC_006.derivFilt   
	false;...	LIC_006.MVinput     
	false;... 	PIC_001.AWBTintterm 
	true;...	PIC_001.derivFilt   
	false;...	PIC_001.MVinput     
	false;... 	TIC_001.AWBTintterm 
	true;...	TIC_001.derivFilt   
	false;...	TIC_001.MVinput     
	false;... 	TIC_002.AWBTintterm 
	true;...	TIC_002.derivFilt   
	false...	TIC_002.MVinput     
%	true ...	Vg,flare
	];
fixed_x	= x(idx==true);
free_x 	= x(idx==false);

%----- Fixed and free inputs:
idu = [...
	true; ...	feed.Qf_in  
	true; ...	feed.Rgo 
	true; ...	feed.BSW_in
	true; ...	feed.T_in  
	true; ...	LV_002.P_out
	true; ...	PV_001.P_out
	true; ...	stk.P
	true; ...	P_001.DP_C
	true; ...	P_001.DP_H
	true; ...	P_002.DP
	true; ...	dslt.Qw
	true; ...	dslt.T
	true; ...	UVRC.cool(1).DP
	true; ...	UVRC.cool(2).DP
	true ...	P_003.DP
	];
fixed_u = u(idu==true);
free_u 	= u(idu==false);

%----- Fixed outputs:
idy 	= zeros(37,1);
%{
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
%}
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
	
	res = uotr_ctrl([],x,u,1,[],par);
	res = res(1:end-1);	%Removing dVg,flaredt from set
	res = res([(1:25)';
		28; ... Removing UVRC.TIC(1) dAWBTIntTermkdt and dderivFiltdt from set
		31; ... Removing UVRC.TIC(2) dAWBTIntTermkdt and dderivFiltdt from set
		34; ... Removing UVRC.MCC dAWBTIntTermkdt and dderivFiltdt from set
		35;36; ...
		39; ... Removing LIC_001 dAWBTIntTermkdt and dderivFiltdt from set
		42; ... Removing LIC_002 dAWBTIntTermkdt and dderivFiltdt from set
		45; ... Removing LIC_003 dAWBTIntTermkdt and dderivFiltdt from set
		48; ... Removing LIC_004 dAWBTIntTermkdt and dderivFiltdt from set
		51; ... Removing LIC_005 dAWBTIntTermkdt and dderivFiltdt from set
		54; ... Removing LIC_006 dAWBTIntTermkdt and dderivFiltdt from set
		57; ... Removing PIC_001 dAWBTIntTermkdt and dderivFiltdt from set
		60; ... Removing TIC_001 dAWBTIntTermkdt and dderivFiltdt from set
		63 	... Removing TIC_002 dAWBTIntTermkdt and dderivFiltdt from set
		]);		
	if ~isempty(yss)
		y = uotr([],x,u,3,[],par);
		res = [res;(yss-y(idy==true))];
	end
end