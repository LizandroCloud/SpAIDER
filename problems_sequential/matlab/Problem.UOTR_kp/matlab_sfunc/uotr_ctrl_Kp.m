function [sys,x0,str,ts,simStateCompliance]=uotr_ctrl_Kp(t,x,u,flag,cinit,par)
%function [sys,x0,str,ts,simStateCompliance]=uotr_ctrl(t,x,ut,u,w,flag,cinit,par)
%uotr_ctrl_sp - Oil treatment unit non-linear model S-Function with PID
%controllers proportional gains as input variables.
%
%This model describes the dynamic behavior of an upstream topside
%facility oil treatment unit, using the black-oil approach with PID
%controllers set-points for process variables such as levels, pressures 
%and temperatures as input variables.
%
%[sys,x0,str,ts,simStateCompliance]=uotr_ctrl_Kp(t,x,u,flag,cinit,par)
%
%sys is the S-function output, which depends on which flag is 
%being inputted by simulink. x0 are the states initial condi-
%tions. See the S-function Level-1 MATLAB reference page for more
%details about str, ts and simStateCompliance S-Function ouputs.
%
%t is the independent variable time. x are the state variables.
%u are the input variables, combined into one vector of manipula-
%ted and disturbs. flag is the Simulink variable inputted to get a
%specific output from the S-function. cinit are variable set for 
%state initial condition computation. par are the model parameters.
%
%----- Model parameters ---------------------------------------------
%black-oil physical properties:
%d60:  					Oil specific gravity              		[]
%Tbn:  					Oil normal boiling point            	[K]
%Cws: 					Water salt concentration            	[g/L]
%dg:  					Gas specific gravity                	[]
%		
%triphase_separator par.SG_001:		
%SG_001.LCS:     		Separation chamber (CS) length      	[m]
%SG_001.LCL:     		Lights chamber (CL) length          	[m]
%SG_001.D:       		Vessel diameter                     	[m]
%SG_001.hweir:   		Weir height               				[m]
%		
%valv_liq par.LV_001:		
%LV_001.Cv:     		Valve constant US units    				[gpm]
%LV_001.FL:     		Liq. pressure recovery factor   		[]
%		
%valv_water par.LV_002:		
%LV_002.Cv:     		Valve constant US units     			[gpm]
%LV_002.FL:     		Liq. pressure recovery factor   		[]
%
%valv_gas par.PV_001:		
%PV_001.Cv:     		Valve constant US units      			[gpm]
%PV_001.Fgm:    		Specific heat cap. ratio factor  		[]
%PV_001.DPrT:   		Pressure drop ratio factor       		[]
%		
%hexch_liq par.P_001:		
%P_001.UA:      		Global heat transfer coefficient    	[kW/K]
%		
%degasser par.V_TO_001:		
%V_TO_001.C:     		Vessel length                     		[m]
%V_TO_001.D:     		Vessel diameter          				[m]
%		
%oil_treater par.TO_001:		
%TO_001.C:     			Vessel length                     		[m]
%TO_001.D:     			Vessel diameter          				[m]
%TO_001.tRL50:   		Res. time for water half-removal   		[s]
%TO_001.gamma:   		Water removal eff. correl. 2nd par		[]
%		
%valv_liq par.LV_003:		
%LV_003.Cv:     		Valve constant US units     			[gpm]
%LV_003.FL:     		Liq. pressure recovery factor   		[]
%		
%valv_water par.LV_004:		
%LV_004.Cv:     		Valve constant US units     			[gpm]
%LV_004.FL:     		Liq. pressure recovery factor   		[]
%		
%pump_water par.B_001:		
%B_001.h0:      		Shut-off head                       	[m]
%B_001.A:       		Pressure-flow pump relation coeff.		[s^2/m^5]

%degasser par.V_TO_002:		
%V_TO_002.C:     		Vessel length                     		[m]
%V_TO_002.D:     		Vessel diameter          				[m]
%		
%oil_treater par.TO_002:		
%TO_002.C:     			Vessel length                     		[m]
%TO_002.D:     			Vessel diameter          				[m]
%TO_002.tRL50:   		Res. time for water half-removal   		[s]
%TO_002.gamma:   		Water removal eff. correl. 2nd par		[]
%		
%valv_liq par.LV_005:		
%LV_005.Cv:     		Valve constant US units     			[gpm]
%LV_005.FL:     		Liq. pressure recovery factor   		[]
%		
%valv_water par.LV_006:		
%LV_006.Cv:     		Valve constant US units     			[gpm]
%LV_006.FL:     		Liq. pressure recovery factor   		[]
%		
%pump_water par.B_002:		
%B_002.h0:      		Shut-off head                       	[m]
%B_002.A:       		Pressure-flow pump relation coeff.		[s^2/m^5]
%		
%pump_liq par.B_003:		
%B_003.h0:      		Shut-off head                       	[m]
%B_003.A:       		Pressure-flow pump relation coeff.		[s^2/m^5]
%	
%vapor recovery unit par.UVRC	
%UVRC.Nr_Max:			Maximum shitf assemble rotation     	[rpm]
%UVRC.nS:				Number of compression stages 			[]
%UVRC.comp(k)			k-th compressor model parameters
%UVRC.ascv(k)			k-th anti-surge control valve parameters
%
%scrubber_gas V_001
%V_001.C:       		Vessel length                            [m]
%V_001.D:       		Vessel diameter                          [m]
%--------------------------------------------------------------------
%----- States -------------------------------------------------------
%x(1):  SG_001.hWCS   				Water phase level in CS          	[m]
%x(2):  SG_001.BSWLCS				Basic sediments & water in CS 		[]
%x(3):  SG_001.hCS      			Total level in CS                  	[m]
%x(4):  SG_001.hCL      			Liquid level in CL                  [m]
%x(5):  SG_001.BSWLCL				Basic sediments & water in CL 		[]
%x(6):  SG_001.P        			Pressure                            [kPa]
%x(7):  SG_001.T        			Temperature                         [K]
%x(8):  V_TO_001.hL     			Liquid level                       	[m]
%x(9):  V_TO_001.BSW				Basic sediments & water 			[]
%x(10): V_TO_001.P      			Pressure                            [kPa]
%x(11): V_TO_001.T      			Temperature                         [K]
%x(12): TO_001.hW      				Water phase level                   [m]
%x(13): TO_001.T       				Temperature                         [K]
%x(14): V_TO_002.hL     			Liquid level                       	[m]
%x(15): V_TO_002.BSW				Basic sediments & water 			[]
%x(16): V_TO_002.P      			Pressure                            [kPa]
%x(17): V_TO_002.T      			Temperature                         [K]
%x(18): TO_002.hW      				Water phase level                   [m]
%x(19): TO_002.T       				Temperature                         [K]
%x(20): UVRC.scrb(1).P 				Pressure 							[kPa]
%x(21): UVRC.scrb(1).T  			Temperature                         [K]
%x(22): UVRC.scrb(2).P 				Pressure 							[kPa]
%x(23): UVRC.scrb(2).T  			Temperature                         [K]
%x(24): UVRC.TIC(1).AWBTintterm 	Integral action with AWBT  			[]
%x(25): UVRC.TIC(2).AWBTintterm 	Integral action with AWBT  			[]
%x(26): UVRC.MCC.AWBTintterm 		Integral action with AWBT  			[]
%x(27): V_001.P 					Pressure 							[kPa]
%x(28): V_001.T  					Temperature                         [K]
%x(29): LIC_001.AWBTintterm  		Integral action with AWBT  			[]
%x(30): LIC_002.AWBTintterm  		Integral action with AWBT  			[]
%x(31): LIC_003.AWBTintterm  		Integral action with AWBT  			[]
%x(32): LIC_004.AWBTintterm  		Integral action with AWBT  			[]
%x(33): LIC_005.AWBTintterm  		Integral action with AWBT  			[]
%x(34): LIC_006.AWBTintterm  		Integral action with AWBT  			[]
%x(35): PIC_001.AWBTintterm  		Integral action with AWBT  			[]
%x(36): TIC_001.AWBTintterm 		Integral action with AWBT  			[]
%x(37): TIC_002.AWBTintterm 		Integral action with AWBT  			[]
%x(38): Vg,flare       				Totalized flaring gas std. vol.     [SMm3]
%--------------------------------------------------------------------
%----- Inputs -------------------------------------------------------
%u(1):  feed.Qf_in     		Liq. std. vol. flow rate    		[Sm^3/d]
%u(2):  feed.Rgo     		Gas-oil ratio     					[Sm^3/Sm^3]
%u(3):  feed.BSW_in    		Basic sediments & water       		[]
%u(4):  feed.T_in      		Temperature                   		[K]
%u(5):  LV_002.P_out    	Outlet pressure            			[kPa]
%u(6):  PV_001.P_out    	Outlet pressure     				[kPa]
%u(7): 	stk.P    			Stock tank pressure            		[kPa]
%u(8): 	P_001.DP_C    		Cold side pressure drop				[kPa]
%u(9): 	P_001.DP_H    		Hot side pressure drop				[kPa]
%u(10): P_002.DP    		Pressure drop						[kPa]
%u(11): dslt.Qw    	 		Desalt water liq. std. vol. flow 	[Sm^3/d]
%u(12): dslt.T    	 		Desalt water temperature			[K]
%u(13):	UVRC.cool(1).DP   	Pressure drop        				[kPa]
%u(14):	UVRC.cool(2).DP   	Pressure drop        				[kPa]
%u(15):	P_003.DP     		Pressure drop 						[kPa]
%u(15):	P_003.DP     		Pressure drop 						[kPa]
%u(16): LIC_001.KC     		Proportional gain                   []
%u(17): LIC_002.KC     		Proportional gain                   []  
%u(18): LIC_003.KC     		Proportional gain                   []
%u(19): LIC_004.KC     		Proportional gain                   []
%u(20): LIC_005.KC     		Proportional gain                   []
%u(21): LIC_006.KC     		Proportional gain                   []
%u(22): PIC_001.KC     		Proportional gain                   []
%u(23): TIC_001.KC     		Proportional gain                   []
%u(24): TIC_002.KC     		Proportional gain                   []
%--------------------------------------------------------------------
%------ Outputs -----------------------------------------------------
%y(1):  SG_001.hWCS   		Water phase level in CS          	[m]
%y(2):  SG_001.BSWLCS   	Basic sediments & water in CS 		[]
%y(3):  SG_001.hCS      	Total level in CS                  	[m]
%y(4):  SG_001.hCL      	Liquid level in CL                  [m]
%y(5):  SG_001.BSWLCL		Basic sediments & water in CL 		[]
%y(6):  SG_001.P        	Pressure                            [kPa]
%y(7):  SG_001.T        	Temperature                         [K]
%y(8):  V_TO_001.hL     	Liquid level                       	[m]
%y(9):  V_TO_001.BSW		Basic sediments & water 			[]
%y(10): V_TO_001.P      	Pressure                            [kPa]
%y(11): V_TO_001.T      	Temperature                         [K]
%y(12): TO_001.hW      		Water phase level                   [m]
%y(13): TO_001.BSW      	Basic sediments & water             []
%y(14): TO_001.T       		Temperature                         [K]
%y(15): V_TO_002.hL     	Liquid level                       	[m]
%y(16): V_TO_002.BSW		Basic sediments & water 			[]
%y(17): V_TO_002.P      	Pressure                            [kPa]
%y(18): V_TO_002.T      	Temperature                         [K]
%y(19): TO_002.hW      		Water phase level                   [m]
%y(20): TO_002.BSW      	Basic sediments & water             []
%y(21): TO_002.T       		Temperature                         [K]
%y(22): V_001.P 			Pressure 							[kPa]
%y(23): V_001.T  			Temperature                         [K]
%y(24): P_002.T       		Temperature                         [K]
%y(25): P_003.T       		Temperature                         [K]
%y(26):	UVRC.scrb(1).T 		Temperature                         [K]
%y(27):	UVRC.scrb(2).T  	Temperature                         [K]
%y(28):	UVRC.comp(1).Qg 	Gas std. vol. flow rate        		[Sm^3/d]
%y(29):	UVRC.comp(2).Qg  	Gas std. vol. flow rate        		[Sm^3/d]
%y(30): stk.Qf    			Crude oil liq. std. vol. flow rate 	[Sm^3/d]
%y(31): stk.Qg   			Crude oil lib. gas std. vol. flow 	[Sm^3/d]
%y(32): stk.BSW   			Crude oil basic sediments & water  	[]
%y(33): stk.T         		Crude oil temperature              	[K]
%y(34): SG_001.Qw    		Water std. vol. flow rate      		[Sm^3/d]
%y(35): TO_001.Qw    		Water std. vol. flow rate      		[Sm^3/d]
%y(36): TO_002.Qw    		Water std. vol. flow rate      		[Sm^3/d]
%y(37): SG_001.Qg    		Gas std. vol. flow rate        		[Sm^3/d]
%y(38): V_TO_001.Qg    		Gas std. vol. flow rate        		[Sm^3/d]
%y(39): V_TO_002.Qg    		Gas std. vol. flow rate        		[Sm^3/d]
%--------------------------------------------------------------------
%------ Initial condition set ---------------------------------------
%cinit = x
%--------------------------------------------------------------------

%2024-06-09. Version 0. Coded by Diego Trica (diegotrica@gmail.com)

%Blackoil object class property package. On the main script the
%user must declare PP as global and PP=blackoil.
%global PP

%----- Parameters:
M_001.nF = 1;
M_001.nW = 2;

M_001.par 		= {par.d60,par.Tbn,par.Cws(1),par.dg(1),M_001.nF,M_001.nW};
SG_001.par 		= [par.d60,par.Tbn,par.Cws(1),par.dg(1),par.SG_001];
LV_001.par 		= [par.d60,par.Tbn,par.Cws(1),par.dg(2),par.LV_001];
LV_002.par 		= [par.Cws(1),par.LV_002];
PV_001.par 		= [par.dg(1),par.PV_001];
FE_001.par 		= [par.dg(1),1e4,par.PV_001(2:3)];
P_001.par 		= {par.d60,par.Tbn,par.Cws(1),par.dg(2),par.P_001};
P_002.par 		= {par.d60,par.Tbn,par.Cws(1),par.dg(2)};
V_TO_001.par 	= [par.d60,par.Tbn,par.Cws(1),par.dg(2),par.V_TO_001];
TO_001.par		= [par.d60,par.Tbn,par.Cws(1),par.dg(2),par.TO_001];
LV_003.par 		= [par.d60,par.Tbn,par.Cws(1),par.dg(3),par.LV_003];
B_001.par 		= [par.Cws(1),par.B_001];
LV_004.par 		= [par.Cws(1),par.LV_004];
FE_002.par 		= {par.dg(2),1e4,0.8,0.7};
D_001.par 		= {par.d60,par.Tbn,par.Cws(2),par.dg(3),1,1};
V_TO_002.par 	= [par.d60,par.Tbn,par.Cws(2),par.dg(3),par.V_TO_002];
TO_002.par		= [par.d60,par.Tbn,par.Cws(2),par.dg(3),par.TO_002];
LV_005.par 		= [par.d60,par.Tbn,par.Cws(2),par.dg(3),par.LV_005];
B_002.par 		= [par.Cws(2),par.B_002];
B_003.par 		= [par.d60,par.Tbn,par.Cws(2),par.dg(3),par.B_003];
LV_006.par 		= [par.Cws(2),par.LV_006];
FE_003.par 		= {par.dg(3),1e4,0.8,0.7};
UVRC.par.d60 	= [par.d60,par.d60];
UVRC.par.dg 	= par.dg(2:3);
UVRC.par.Nr_Max = par.UVRC.Nr_Max;
UVRC.par.nS		= par.UVRC.nS;
UVRC.par.Qout_Max = par.UVRC.Qout_Max;
UVRC.par.scrb	= par.UVRC.scrb;
UVRC.par.comp	= par.UVRC.comp;
UVRC.par.ascv	= par.UVRC.ascv;
G_001.par		= {par.dg(1),2};
P_003.par 		= {par.dg(1)};
V_001.par		= [par.d60,par.dg(1),par.V_001];
LIT_001.par		= {0.00,3.00}; 	%{RangeMin,RangeMax} [m]
LIC_001.par 	= {'Direct'};	%{action}
LIT_002.par		= {0.00,3.00}; 	%{RangeMin,RangeMax} [m]
LIC_002.par 	= {'Direct'};	%{action}
LIT_003.par		= {0.00,3.15}; 	%{RangeMin,RangeMax} [m]
LIC_003.par 	= {'Direct'};	%{action}
LIT_004.par		= {0.00,4.30}; 	%{RangeMin,RangeMax} [m]
LIC_004.par 	= {'Direct'};	%{action}
LIT_005.par		= {0.00,3.15}; 	%{RangeMin,RangeMax} [m]
LIC_005.par 	= {'Direct'};	%{action}
LIT_006.par		= {0.00,4.30}; 	%{RangeMin,RangeMax} [m]
LIC_006.par 	= {'Direct'};	%{action}
PIT_001.par		= {0.00,3e3}; 	%{RangeMin,RangeMax} [kPa]
PIC_001.par 	= {'Direct'};	%{action}
TIT_001.par		= {270,400}; 	%{RangeMin,RangeMax} [K]
TIC_001.par 	= {'Reverse'};	%{action}
TIT_002.par		= {270,400}; 	%{RangeMin,RangeMax} [K]
TIC_002.par 	= {'Direct'};	%{action}

%----- States:
%triphase_separator SG_001
SG_001.x 			= x(1:7);
SG_001.hWCS  		= SG_001.x(1);
SG_001.BSWLCS 		= SG_001.x(2);
SG_001.hCS  		= SG_001.x(3);
SG_001.hCL  		= SG_001.x(4);
SG_001.BSWLCL 		= SG_001.x(5);
SG_001.P       		= SG_001.x(6);
SG_001.T			= SG_001.x(7);

%degasser V_TO_001
V_TO_001.x			= x(8:11);
V_TO_001.hL 		= V_TO_001.x(1);
V_TO_001.BSW		= V_TO_001.x(2);
V_TO_001.P  		= V_TO_001.x(3);
V_TO_001.T  		= V_TO_001.x(4);

%oil_treater TO_001
TO_001.x			= x(12:13);
TO_001.hW 			= TO_001.x(1);
TO_001.T  			= TO_001.x(2);

%degasser V_TO_002
V_TO_002.x			= x(14:17);
V_TO_002.hL 		= V_TO_002.x(1);
V_TO_002.BSW		= V_TO_002.x(2);
V_TO_002.P  		= V_TO_002.x(3);
V_TO_002.T  		= V_TO_002.x(4);

%oil_treater TO_002
TO_002.x			= x(18:19);
TO_002.hW 			= TO_002.x(1);
TO_002.T  			= TO_002.x(2);

%vapor recovery unit UVRC
UVRC.x					= x(20:26);
UVRC.scrb(1).P 			= UVRC.x(1);
UVRC.scrb(1).T 			= UVRC.x(2);
UVRC.scrb(2).P 			= UVRC.x(3);
UVRC.scrb(2).T 			= UVRC.x(4);
UVRC.TIC(1).AWBTintterm = UVRC.x(5);
UVRC.TIC(2).AWBTintterm = UVRC.x(6);
UVRC.MCC.AWBTintterm 	= UVRC.x(7);

%scrubber_gas V_001
V_001.x				= x(27:28);
V_001.P 			= V_001.x(1);
V_001.T 			= V_001.x(2);

%PI_pos LIC_001
LIC_001.x			= x(29);
LIC_001.AWBTintterm = LIC_001.x(1);

%PI_pos LIC_002
LIC_002.x			= x(30);
LIC_002.AWBTintterm = LIC_002.x(1);

%PI_pos LIC_003
LIC_003.x			= x(31);
LIC_003.AWBTintterm = LIC_003.x(1);

%PI_pos LIC_004
LIC_004.x			= x(32);
LIC_004.AWBTintterm = LIC_004.x(1);

%PI_pos LIC_005
LIC_005.x			= x(33);
LIC_005.AWBTintterm = LIC_005.x(1);

%PI_pos LIC_006
LIC_006.x			= x(34);
LIC_006.AWBTintterm = LIC_006.x(1);

%PI_pos PIC_001
PIC_001.x			= x(35);
PIC_001.AWBTintterm = PIC_001.x(1);

%PI_pos TIC_001
TIC_001.x			= x(36);
TIC_001.AWBTintterm = TIC_001.x(1);

%PI_pos TIC_002
TIC_002.x			= x(37);
TIC_002.AWBTintterm = TIC_002.x(1);

%----- Inputs:
feed.Qf_in			= u(1);
feed.Rgo			= u(2);
feed.BSW_in			= u(3);
feed.T_in			= u(4);
LV_002.P_out		= u(5);
PV_001.P_out		= u(6);
stk.P  				= u(7);
P_001.DP_C  		= u(8);
P_001.DP_H  		= u(9);
P_002.DP  			= u(10);
dslt.Qw    			= u(11);
dslt.T     			= u(12);
UVRC.cool(1).DP   	= u(13);
UVRC.cool(2).DP   	= u(14);
P_003.DP			= u(15);
% LIC_001.KC 			= u(16);
% LIC_002.KC 			= u(17);
% LIC_003.KC 			= u(18);
% LIC_004.KC 			= u(19);
% LIC_005.KC 			= u(20);
% LIC_006.KC 			= u(21);
% PIC_001.KC 			= u(22);
% TIC_001.KC 			= u(23);
% TIC_002.KC 			= u(24);

LIC_001.KC 			= u(16);
LIC_003.KC 			= u(17);
LIC_005.KC 			= u(18);

%----- Control parameters:
%LIC_001: SG_001.hCL liquid level control
LIC_001.mode		= 0;	%auto
LIT_001.SPinput 	= 1.50;	%[m]
% LIC_001.KC			= 3.0;%4.0;
LIC_001.TI			= 900.0; %600.0;%[s]
LIC_001.TD			= 0.0;	%[s] (PI controller)

%LIC_002: SG_001.hWCS liquid level control
LIC_002.mode		= 0;	%auto
LIT_002.SPinput 	= 0.15;	%[m]
LIC_002.KC			= 3.0;%4.0; 
LIC_002.TI			= 900.0; %600.0;%[s]
LIC_002.TD			= 0.0;	%[s] (PI controller)

%LIC_003: V_TO_001.hL liquid level control
LIC_003.mode		= 0;	%auto
LIT_003.SPinput 	= 1.00;	%[m]
% LIC_003.KC			= 3.0;%4.0;
LIC_003.TI			= 900.0;%[s]600.0;%[s]
LIC_003.TD			= 0.0;	%[s] (PI controller)

%LIC_003: TO_001.hW liquid level control
LIC_004.mode		= 0;	%auto
LIT_004.SPinput 	= 2.00;	%[m]
LIC_004.KC			= 3.0;%4.0;
LIC_004.TI			= 900.0;%[s]600.0;%[s]
LIC_004.TD			= 0.0;	%[s] (PI controller)

%LIC_005: V_TO_002.hL liquid level control
LIC_005.mode		= 0;	%auto
LIT_005.SPinput 	= 1.00;	%[m]
% LIC_005.KC			= 3.0; %4.0;
LIC_005.TI			= 900.0;%[s]600.0;%[s]
LIC_005.TD			= 0.0;	%[s] (PI controller)

%LIC_006: TO_002.hW liquid level control
LIC_006.mode		= 0;	%auto
LIT_006.SPinput 	= 2.00;	%[m]
LIC_006.KC			= 3.0;%3.0; %4.0;
LIC_006.TI			= 900;%900.0;%[s]600.0;%[s]
LIC_006.TD			= 0.0;	%[s] (PI controller)

%PIC_001: V_001.P pressure control
PIC_001.mode		= 0;	%auto
PIT_001.SPinput 	= 2101.325;	%[kPa]
PIC_001.KC			= 3.0;%2.5;%1.0;
PIC_001.TI			= 60;%45.0;%60.0; %[s]
PIC_001.TD			= 0.0;	%[s] (PI controller)

%TIC_001: V_TO_001.T temperature control
TIC_001.mode		= 0;	%auto
TIT_001.SPinput 	= 360;	%[K]
TIC_001.KC			= 1.0;%1.0;
TIC_001.TI			= 120;%60.0; %[s]
TIC_001.TD			= 0.0;	%[s] (PI controller)

%TIC_002: V_001.T temperature control
TIC_002.mode		= 0;	%auto
TIT_002.SPinput 	= 310;	%[K]
TIC_002.KC			= 1.0;%1.0;
TIC_002.TI			= 120;%60.0; %[s]
TIC_002.TD			= 0.0;	%[s] (PI controller)

%----- Internal calculations:
%feed stream Gas std vol. flow rate [Sm^3/d]
feed.Qg_in = feed.Qf_in*(1-feed.BSW_in)*feed.Rgo;	

%transmitter LIT_001
LIT_001.PVinput = SG_001.hCL;
LIT_001.u = [...
	LIT_001.PVinput; ...
	LIT_001.SPinput];
LIT_001.y = transmitter(t,[],LIT_001.u,3,[],LIT_001.par);
LIT_001.PVoutput = LIT_001.y(1);
LIT_001.SPoutput = LIT_001.y(2);

%PI_pos LIC_001
LIC_001.PV 	= LIT_001.PVoutput;
LIC_001.SP	= LIT_001.SPoutput;
LIC_001.u = [...
	LIC_001.mode; ...
	LIC_001.PV; ...
	LIC_001.SP; ...
	LIC_001.KC; LIC_001.TI; LIC_001.TD];
LIC_001.y = PI_pos(t,LIC_001.x,LIC_001.u,3,[],LIC_001.par);
LIC_001.MVoutput = LIC_001.y(1);

%valv_liq LV_001:
LV_001.P_in 	= SG_001.P;
LV_001.P_out	= V_TO_001.P + P_002.DP + P_001.DP_C;
LV_001.c 		= LIC_001.MVoutput;
LV_001.u = [...
    SG_001.BSWLCL; ...
    LV_001.P_in; ...
    SG_001.T; ...
    LV_001.P_out; ...
    LV_001.c];
LV_001.y = valv_liq(t,[],LV_001.u,3,[],LV_001.par);
SG_001.Qf 		= LV_001.y(1);
LV_001.Qg 		= LV_001.y(2);
LV_001.BSW 		= LV_001.y(3);
LV_001.T 		= LV_001.y(4);

%transmitter LIT_002
LIT_002.PVinput = SG_001.hWCS;
LIT_002.u = [...
	LIT_002.PVinput; ...
	LIT_002.SPinput];
LIT_002.y = transmitter(t,[],LIT_002.u,3,[],LIT_002.par);
LIT_002.PVoutput = LIT_002.y(1);
LIT_002.SPoutput = LIT_002.y(2);

%PI_pos LIC_002
LIC_002.PV 	= LIT_002.PVoutput;
LIC_002.SP	= LIT_002.SPoutput;
LIC_002.u = [...
	LIC_002.mode; ...
	LIC_002.PV; ...
	LIC_002.SP; ...
	LIC_002.KC; LIC_002.TI; LIC_002.TD];
LIC_002.y = PI_pos(t,LIC_002.x,LIC_002.u,3,[],LIC_002.par);
LIC_002.MVoutput = LIC_002.y(1);

%valv_water LV_002:
LV_002.P_in = SG_001.P;
LV_002.c 	= LIC_002.MVoutput;
LV_002.u = [...
    LV_002.P_in; ...
    SG_001.T; ...
    LV_002.P_out; ...
    LV_002.c];
LV_002.y = valv_water(t,[],LV_002.u,3,[],LV_002.par);
SG_001.Qw       = LV_002.y(1);

%Iterative calculation for z=[P_001.Qg_C; P_001.T_outC]
jter=0; maxjter=100; rtol=1e-5;
P_001.Qg_C		= ...
	LV_001.Qg + 1*SG_001.Qf*(1-LV_001.BSW);	%initial guess for P_001.Qg_C
P_001.T_outC	= LV_001.T + 20;			%initial guess for P_001.T_outC
while jter<maxjter
    jter=jter+1;
    z0 = [P_001.Qg_C; P_001.T_outC];

%transmitter TIT_001
TIT_001.PVinput = V_TO_001.T;
TIT_001.u = [...
	TIT_001.PVinput; ...
	TIT_001.SPinput];
TIT_001.y = transmitter(t,[],TIT_001.u,3,[],TIT_001.par);
TIT_001.PVoutput = TIT_001.y(1);
TIT_001.SPoutput = TIT_001.y(2);

%PI_pos TIC_001
TIC_001.PV 	= TIT_001.PVoutput;
TIC_001.SP	= TIT_001.SPoutput;
TIC_001.u = [...
	TIC_001.mode; ...
	TIC_001.PV; ...
	TIC_001.SP; ...
	TIC_001.KC; TIC_001.TI; TIC_001.TD];
TIC_001.y = PI_pos(t,TIC_001.x,TIC_001.u,3,[],TIC_001.par);
TIC_001.MVoutput = TIC_001.y(1);

%Production heater P_002
P_002.Qf_in 	= SG_001.Qf;
P_002.Qg_in 	= P_001.Qg_C;
P_002.BSW_in 	= LV_001.BSW;
P_002.P_in 		= V_TO_001.P + P_002.DP;
P_002.T_in 		= P_001.T_outC;
P_002.P_out 	= V_TO_001.P;
%P_002.Q_in 		= (5 + (25 - 5)*TIC_001.MVoutput)*1e3; %[kW]
P_002.Q_in 		= (0 + (25 - 0)*TIC_001.MVoutput)*1e3; %[kW]
P_002.u = [...
    P_002.Qf_in; ...
	P_002.Qg_in; ...
	P_002.BSW_in; ...
    P_002.P_in; ...
    P_002.T_in; ...
    P_002.P_out; ...
    P_002.Q_in];
P_002.y = heater_liq(t,[],P_002.u,3,[],P_002.par);
P_002.Qf        = P_002.y(1);
P_002.Qg        = P_002.y(2);
P_002.BSW       = P_002.y(3);
P_002.T         = P_002.y(4);

%Iterative calculation for y=[TO_001.Qf; TO_001.Qw]
TO_001.BSW_in 	= V_TO_001.BSW;
TO_001.P_in 	= V_TO_001.P; 
TO_001.T_in 	= V_TO_001.T;

iter=0; maxiter=100;
if jter == 1
    TO_001.Qf		= P_002.Qf*(1-V_TO_001.BSW);	%initial guess for TO_001.Qf
    TO_001.Qw		= P_002.Qf*V_TO_001.BSW;		%initial guess for TO_001.Qw
end
while iter<maxiter
    iter=iter+1;
    y0 = [TO_001.Qf; TO_001.Qw];
	
	%oil_treater TO_001:
	TO_001.u = [...
        TO_001.BSW_in; ...
        TO_001.P_in; ...
        TO_001.T_in; ...
        TO_001.Qf; ...
        TO_001.Qw];
	TO_001.y = oil_treater(t,TO_001.x,TO_001.u,3,[],TO_001.par);
	TO_001.BSW      = TO_001.y(2);
	TO_001.Qf_in    = TO_001.y(3);
	TO_001.T        = TO_001.y(4);	
	
	%transmitter LIT_003
	LIT_003.PVinput = V_TO_001.hL;
	LIT_003.u = [...
		LIT_003.PVinput; ...
		LIT_003.SPinput];
	LIT_003.y = transmitter(t,[],LIT_003.u,3,[],LIT_003.par);
	LIT_003.PVoutput = LIT_003.y(1);
	LIT_003.SPoutput = LIT_003.y(2);

	%PI_pos LIC_003
	LIC_003.PV 	= LIT_003.PVoutput;
	LIC_003.SP	= LIT_003.SPoutput;
	LIC_003.u = [...
		LIC_003.mode; ...
		LIC_003.PV; ...
		LIC_003.SP; ...
		LIC_003.KC; LIC_003.TI; LIC_003.TD];
	LIC_003.y = PI_pos(t,LIC_003.x,LIC_003.u,3,[],LIC_003.par);
	LIC_003.MVoutput = LIC_003.y(1);	
	
	%valv_liq LV_003:
	LV_003.BSW_in   = TO_001.BSW;
	LV_003.c 		= LIC_003.MVoutput;
	LV_003.u = [...
        LV_003.BSW_in; ...
        V_TO_001.P; ...
        TO_001.T; ...
        V_TO_002.P; ...
        LV_003.c];
	LV_003.y = valv_liq(t,[],LV_003.u,3,[],LV_003.par);
	TO_001.Qf       = LV_003.y(1);
	
	%pump_water B_004:
	B_001.Qw_in     = TO_001.Qw;
	B_001.P_in      = V_TO_001.P;
	B_001.T_in      = TO_001.T;
	B_001.u = [...
        B_001.Qw_in; ...
        B_001.P_in; ...
        B_001.T_in];
	B_001.y = pump_water(t,[],B_001.u,3,[],B_001.par);
	B_001.P_out     = B_001.y(1);

	%transmitter LIT_004
	LIT_004.PVinput = TO_001.hW;
	LIT_004.u = [...
		LIT_004.PVinput; ...
		LIT_004.SPinput];
	LIT_004.y = transmitter(t,[],LIT_004.u,3,[],LIT_004.par);
	LIT_004.PVoutput = LIT_004.y(1);
	LIT_004.SPoutput = LIT_004.y(2);

	%PI_pos LIC_004
	LIC_004.PV 	= LIT_004.PVoutput;
	LIC_004.SP	= LIT_004.SPoutput;
	LIC_004.u = [...
		LIC_004.mode; ...
		LIC_004.PV; ...
		LIC_004.SP; ...
		LIC_004.KC; LIC_004.TI; LIC_004.TD];
	LIC_004.y = PI_pos(t,LIC_004.x,LIC_004.u,3,[],LIC_004.par);
	LIC_004.MVoutput = LIC_004.y(1);

	%valv_water LV_004:
	LV_004.P_in 	= B_001.P_out;
	LV_004.T_in 	= TO_001.T;
	LV_004.P_out 	= SG_001.P;
	LV_004.c 		= LIC_004.MVoutput;
	LV_004.u = [...
        LV_004.P_in; ...
        LV_004.T_in; ...
        LV_004.P_out; ...
        LV_004.c];
	LV_004.y = valv_water(t,[],LV_004.u,3,[],LV_004.par);
	TO_001.Qw       = LV_004.y(1);

	y = [TO_001.Qf; TO_001.Qw];
	rdev=abs((y-y0)./y0);
    if rdev<rtol
        break
    end
end

%mixer_W3ph D_001
D_001.F 	= struct('Qf',LV_003.y(1),'Qg',LV_003.y(2),'BSW',LV_003.y(3),'T',LV_003.y(4));
D_001.W		= struct('Qw',dslt.Qw,'T',dslt.T);
D_001.u = [...
    struct2cell(D_001.F); ...
    struct2cell(D_001.W); ...
    V_TO_002.P];
D_001.u = cell2mat(D_001.u);
D_001.y = mixer_W3ph(t,[],D_001.u,3,[],D_001.par);
D_001.Qf        = D_001.y(1);
D_001.Qg        = D_001.y(2);
D_001.BSW       = D_001.y(3);
D_001.T         = D_001.y(4);

%Iterative calculation for y=[TO_002.Qf; TO_002.Qw]
TO_002.BSW_in 	= V_TO_002.BSW;
TO_002.P_in 	= V_TO_002.P; 
TO_002.T_in 	= V_TO_002.T;

iter=0; maxiter=100;
if jter == 1
    TO_002.Qf		= D_001.Qf*(1-V_TO_002.BSW);	%initial guess for TO_002.Qf
    TO_002.Qw		= D_001.Qf*V_TO_002.BSW;		%initial guess for TO_002.Qw
end
while iter<maxiter
    iter=iter+1;
    y0 = [TO_002.Qf; TO_002.Qw];
	
	%oil_treater TO_002:
	TO_002.u = [...
        TO_002.BSW_in; ...
        TO_002.P_in; ...
        TO_002.T_in; ...
        TO_002.Qf; ...
        TO_002.Qw];
	TO_002.y = oil_treater(t,TO_002.x,TO_002.u,3,[],TO_002.par);
	TO_002.BSW      = TO_002.y(2);
	TO_002.Qf_in    = TO_002.y(3);
	TO_002.T        = TO_002.y(4);	

	%pump_liq B_003:
	B_003.Qf_in     = TO_002.Qf;
	B_003.BSW_in    = TO_002.BSW;
	B_003.P_in      = V_TO_002.P;
	B_003.T_in      = TO_002.T;
	B_003.u = [...
        B_003.Qf_in; ...
		B_003.BSW_in; ...
        B_003.P_in; ...
        B_003.T_in];
	B_003.y = pump_liq(t,[],B_003.u,3,[],B_003.par);
	B_003.P_out     = B_003.y(1);
	
	%transmitter LIT_005
	LIT_005.PVinput = V_TO_002.hL;
	LIT_005.u = [...
		LIT_005.PVinput; ...
		LIT_005.SPinput];
	LIT_005.y = transmitter(t,[],LIT_005.u,3,[],LIT_005.par);
	LIT_005.PVoutput = LIT_005.y(1);
	LIT_005.SPoutput = LIT_005.y(2);

	%PI_pos LIC_005
	LIC_005.PV 	= LIT_005.PVoutput;
	LIC_005.SP	= LIT_005.SPoutput;
	LIC_005.u = [...
		LIC_005.mode; ...
		LIC_005.PV; ...
		LIC_005.SP; ...
		LIC_005.KC; LIC_005.TI; LIC_005.TD];
	LIC_005.y = PI_pos(t,LIC_005.x,LIC_005.u,3,[],LIC_005.par);
	LIC_005.MVoutput = LIC_005.y(1);	
	
	%valv_liq LV_005:
	LV_005.BSW_in   = TO_002.BSW;
	LV_005.P_in 	= B_003.P_out;
	LV_005.T_in 	= TO_002.T;
	LV_005.P_out 	= stk.P + P_001.DP_H;
	LV_005.c		= LIC_005.MVoutput;
	LV_005.u = [...
        LV_005.BSW_in; ...
        LV_005.P_in; ...
        LV_005.T_in; ...
        LV_005.P_out; ...
        LV_005.c];
	LV_005.y = valv_liq(t,[],LV_005.u,3,[],LV_005.par);
	TO_002.Qf       = LV_005.y(1);
	LV_005.Qg 		= LV_005.y(2);
	LV_005.T 		= LV_005.y(4);

	%pump_water B_002:
	B_002.Qw_in     = TO_002.Qw;
	B_002.P_in      = V_TO_002.P;
	B_002.T_in      = TO_002.T;
	B_002.u = [...
        B_002.Qw_in; ...
        B_002.P_in; ...
        B_002.T_in];
	B_002.y = pump_water(t,[],B_002.u,3,[],B_002.par);
	B_002.P_out     = B_002.y(1);

	%transmitter LIT_006
	LIT_006.PVinput = TO_002.hW;
	LIT_006.u = [...
		LIT_006.PVinput; ...
		LIT_006.SPinput];
	LIT_006.y = transmitter(t,[],LIT_006.u,3,[],LIT_006.par);
	LIT_006.PVoutput = LIT_006.y(1);
	LIT_006.SPoutput = LIT_006.y(2);

	%PI_pos LIC_006
	LIC_006.PV 	= LIT_006.PVoutput;
	LIC_006.SP	= LIT_006.SPoutput;
	LIC_006.u = [...
		LIC_006.mode; ...
		LIC_006.PV; ...
		LIC_006.SP; ...
		LIC_006.KC; LIC_006.TI; LIC_006.TD];
	LIC_006.y = PI_pos(t,LIC_006.x,LIC_006.u,3,[],LIC_006.par);
	LIC_006.MVoutput = LIC_006.y(1);

	%valv_water LV_006:
	LV_006.P_in 	= B_002.P_out;
	LV_006.T_in 	= TO_002.T;
	LV_006.P_out 	= SG_001.P;
	LV_006.c		= LIC_006.MVoutput;
	LV_006.u = [...
        LV_006.P_in; ...
        LV_006.T_in; ...
        LV_006.P_out; ...
        LV_006.c];
	LV_006.y = valv_water(t,[],LV_006.u,3,[],LV_006.par);
	TO_002.Qw       = LV_006.y(1);

	y = [TO_002.Qf; TO_002.Qw];
	rdev=abs((y-y0)./y0);
    if rdev<rtol
        break
    end
end

%Oil-oil heat exchanger P_001
P_001.Qf_inC 	= SG_001.Qf;
P_001.Qg_inC 	= LV_001.Qg;
P_001.BSW_inC 	= LV_001.BSW;
P_001.P_inC 	= V_TO_001.P + P_002.DP + P_001.DP_C;
P_001.T_inC 	= LV_001.T;
P_001.P_outC 	= V_TO_001.P + P_002.DP;
P_001.Qf_inH 	= TO_002.Qf;
P_001.Qg_inH 	= LV_005.Qg;
P_001.BSW_inH 	= TO_002.BSW;
P_001.P_inH 	= stk.P + P_001.DP_H;
P_001.P_outH 	= stk.P;
P_001.T_inH 	= LV_005.T;
P_001.u = [...
    P_001.Qf_inC; ...
	P_001.Qg_inC; ...
	P_001.BSW_inC; ...
    P_001.P_inC; ...
    P_001.T_inC; ...
    P_001.P_outC; ...
    P_001.Qf_inH; ...
    P_001.Qg_inH; ...
    P_001.BSW_inH; ...
    P_001.P_inH; ...
    P_001.T_inH; ...
    P_001.P_outH];
P_001.y = hexch_liq(t,[],P_001.u,3,[],P_001.par);
P_001.Qg_C		= P_001.y(2);
P_001.T_outC	= P_001.y(4);

%Crude oil to stock tank stream variables
stk.Qf 		= P_001.y(5);
stk.Qg 		= P_001.y(6);
stk.BSW 	= P_001.y(7);
stk.T 		= P_001.y(8);

	z = [P_001.Qg_C; P_001.T_outC];
	rdev=abs((z-z0)./z0);
    if rdev<rtol
        break
    end
end

%mixer_W3ph M_001
%M_001.F 	= feed;
M_001.F 	= struct('Qf',feed.Qf_in,'Qg',feed.Qg_in,'BSW',feed.BSW_in,'T',feed.T_in);
M_001.W(1)	= struct('Qw',LV_004.y(1),'T',LV_004.y(2));
M_001.W(2)	= struct('Qw',LV_006.y(1),'T',LV_006.y(2));
M_001.u = [...
    struct2cell(M_001.F); ...
    struct2cell(M_001.W(1)); ...
	struct2cell(M_001.W(2)); ...
    SG_001.P];
M_001.u = cell2mat(M_001.u);
M_001.y = mixer_W3ph(t,[],M_001.u,3,[],M_001.par);
SG_001.Qf_in 	= M_001.y(1);
SG_001.Qg_in 	= M_001.y(2);
SG_001.BSW_in 	= M_001.y(3);
SG_001.T_in 	= M_001.y(4);

%valv_gas FE_002:
FE_002.u = [V_TO_001.P; V_TO_001.T; (UVRC.scrb(2).P + UVRC.cool(2).DP); 1.00];
FE_002.y = valv_gas(t,[],FE_002.u,3,[],FE_002.par);
FE_002.Qg		= FE_002.y(1);
FE_002.T		= FE_002.y(2);

%valv_gas FE_003:
FE_003.u = [V_TO_002.P; V_TO_002.T; (UVRC.scrb(1).P + UVRC.cool(1).DP); 1.00];
FE_003.y = valv_gas(t,[],FE_003.u,3,[],FE_003.par);
FE_003.Qg		= FE_003.y(1);
FE_003.T		= FE_003.y(2);

%degasser V_TO_001 Gas std. vol. flow rate [Sm^3/d]
V_TO_001.Qg = FE_002.Qg;

%degasser V_TO_002 Gas std. vol. flow rate [Sm^3/d]
V_TO_002.Qg = FE_003.Qg;

%vapor recovery unit UVRC
UVRC.P_out = (V_001.P + P_003.DP);
UVRC.feed(1).Qg = V_TO_002.Qg;
UVRC.feed(1).T = FE_002.T;
UVRC.feed(2).Qg = V_TO_001.Qg;
UVRC.feed(2).T = FE_003.T;

for k=1:UVRC.par.nS
UVRC.u(2*(k-1)+1)					= UVRC.feed(k).Qg;
UVRC.u(2*(k-1)+2)					= UVRC.feed(k).T;
UVRC.u(2*UVRC.par.nS+(k-1)+1)    	= UVRC.cool(k).DP;
end
UVRC.u(2*UVRC.par.nS+UVRC.par.nS+1)	= UVRC.P_out;
UVRC.u = UVRC.u(:);
UVRC.y = uvrc_ctrl(t,UVRC.x,UVRC.u,3,[],UVRC.par);

for k=1:UVRC.par.nS
UVRC.feed(k).P 		= UVRC.y((k-1)+1);
UVRC.comp(k).Qg		= UVRC.y(UVRC.par.nS+3*(k-1)+1);
UVRC.comp(k).T		= UVRC.y(UVRC.par.nS+3*(k-1)+2);
UVRC.comp(k).W		= UVRC.y(UVRC.par.nS+3*(k-1)+3);
UVRC.stg(k).Qg_th	= UVRC.y(UVRC.par.nS+3*UVRC.par.nS+(k-1)+1);
UVRC.scrb(k).P		= UVRC.y(UVRC.par.nS+3*UVRC.par.nS+2*UVRC.par.nS+2*(k-1)+1);
UVRC.scrb(k).T		= UVRC.y(UVRC.par.nS+3*UVRC.par.nS+2*UVRC.par.nS+2*(k-1)+2);
end

%valv_gas FE_001:
FE_001.P_out = V_001.P + P_003.DP;
FE_001.u = [SG_001.P; SG_001.T; FE_001.P_out; 1.00];
FE_001.y = valv_gas(t,[],FE_001.u,3,[],FE_001.par);
FE_001.Qg		= FE_001.y(1);
FE_001.T		= FE_001.y(2);

%triphase_separator SG_001 Gas std. vol. flow rate [Sm^3/d]
SG_001.Qg = FE_001.Qg;

%mixer_G G_001
G_001.G(1).Qg 		= FE_001.Qg;
G_001.G(2).Qg 		= UVRC.stg(2).Qg_th;
	
G_001.G(1).T 		= FE_001.T;
G_001.G(2).T 		= UVRC.comp(2).T;

G_001.u = [];
for j=1:G_001.par{2}
    G_001.u(2*(j-1)+1) = G_001.G(j).Qg;
    G_001.u(2*(j-1)+2) = G_001.G(j).T;
end
G_001.u = [G_001.u, SG_001.P];
G_001.u = G_001.u(:);
G_001.y = mixer_G(t,[],G_001.u,3,[],G_001.par);
G_001.Qg 		= G_001.y(1);
G_001.T 		= G_001.y(2);

%transmitter TIT_002
TIT_002.PVinput = V_001.T;
TIT_002.u = [...
	TIT_002.PVinput; ...
	TIT_002.SPinput];
TIT_002.y = transmitter(t,[],TIT_002.u,3,[],TIT_002.par);
TIT_002.PVoutput = TIT_002.y(1);
TIT_002.SPoutput = TIT_002.y(2);

%PI_pos TIC_002
TIC_002.PV 	= TIT_002.PVoutput;
TIC_002.SP	= TIT_002.SPoutput;
TIC_002.u = [...
	TIC_002.mode; ...
	TIC_002.PV; ...
	TIC_002.SP; ...
	TIC_002.KC; TIC_002.TI; TIC_002.TD];
TIC_002.y = PI_pos(t,TIC_002.x,TIC_002.u,3,[],TIC_002.par);
TIC_002.MVoutput = TIC_002.y(1);

%cooler_gas P_003
P_003.Qg_in 	= G_001.Qg;
P_003.P_in 		= V_001.P + P_003.DP;
P_003.T_in 		= G_001.T;
P_003.P_out 	= V_001.P;
P_003.Q_out 	= 2.5e3*TIC_002.MVoutput; 	%[kW]
P_003.u = [...
	P_003.Qg_in; ...
	P_003.P_in; ...
	P_003.T_in; ...
	P_003.P_out; ...
	P_003.Q_out];
P_003.y = cooler_gas(t,[],P_003.u,3,[],P_003.par);
P_003.Qg 		= P_003.y(1);
P_003.T 		= P_003.y(2);

%transmitter PIT_001
PIT_001.PVinput = V_001.P;
PIT_001.u = [...
	PIT_001.PVinput; ...
	PIT_001.SPinput];
PIT_001.y = transmitter(t,[],PIT_001.u,3,[],PIT_001.par);
PIT_001.PVoutput = PIT_001.y(1);
PIT_001.SPoutput = PIT_001.y(2);

%PI_pos PIC_001
PIC_001.PV 	= PIT_001.PVoutput;
PIC_001.SP	= PIT_001.SPoutput;
PIC_001.u = [...
	PIC_001.mode; ...
	PIC_001.PV; ...
	PIC_001.SP; ...
	PIC_001.KC; PIC_001.TI; PIC_001.TD];
PIC_001.y = PI_pos(t,PIC_001.x,PIC_001.u,3,[],PIC_001.par);
PIC_001.MVoutput = PIC_001.y(1);

%valv_gas PV_001:	
PV_001.P_in 	= V_001.P;
PV_001.T_in 	= V_001.T;
PV_001.c		= PIC_001.MVoutput;
PV_001.u = [...
    PV_001.P_in; ...
    PV_001.T_in; ...
    PV_001.P_out; ...
    PV_001.c];
PV_001.y = valv_gas(t,[],PV_001.u,3,[],PV_001.par);
PV_001.Qg   = PV_001.y(1);

switch flag,
%----- S-FUNCTION BUILDING PARAMETERS AND INITIAL CONDITIONS --------
    case 0,        
sys(1)=38; %Number of continuous states (P_pos in uvrc)
sys(2)=0;%Number of discrete states
sys(3)=39;%Number of outputs
sys(4)=18;%Number of inputs
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
%----- State derivatives:
%triphase_separator SG_001
SG_001.u = [...
    SG_001.Qf_in; ...
    SG_001.Qg_in; ...
    SG_001.BSW_in; ...
    SG_001.T_in; ...
    SG_001.Qf; ...
    SG_001.Qw; ...
    SG_001.Qg];
SG_001.dxdt = triphase_separator(t,SG_001.x,SG_001.u,1,[],SG_001.par);

%degasser V_TO_001
V_TO_001.Qf_in  = P_002.Qf;
V_TO_001.Qg_in	= P_002.Qg;
V_TO_001.BSW_in	= P_002.BSW;
V_TO_001.T_in	= P_002.T;
V_TO_001.Qf		= TO_001.Qf_in;
V_TO_001.u = [...
    V_TO_001.Qf_in; ...
    V_TO_001.Qg_in; ...
    V_TO_001.BSW_in; ...
    V_TO_001.T_in; ...
    V_TO_001.Qf; ...
    V_TO_001.Qg];
V_TO_001.dxdt = degasser(t,V_TO_001.x,V_TO_001.u,1,[],V_TO_001.par);

%oil_treater TO_001
TO_001.dxdt = oil_treater(t,TO_001.x,TO_001.u,1,[],TO_001.par);

%degasser V_TO_002
V_TO_002.Qf_in  = D_001.Qf;
V_TO_002.Qg_in	= D_001.Qg;
V_TO_002.BSW_in	= D_001.BSW;
V_TO_002.T_in	= D_001.T;
V_TO_002.Qf		= TO_002.Qf_in;
V_TO_002.u = [...
    V_TO_002.Qf_in; ...
    V_TO_002.Qg_in; ...
    V_TO_002.BSW_in; ...
    V_TO_002.T_in; ...
    V_TO_002.Qf; ...
    V_TO_002.Qg];
V_TO_002.dxdt = degasser(t,V_TO_002.x,V_TO_002.u,1,[],V_TO_002.par);

%oil_treater TO_002
TO_002.dxdt = oil_treater(t,TO_002.x,TO_002.u,1,[],TO_002.par);

%vapor recovery unit UVRC_ctrl
UVRC.dxdt = uvrc_ctrl(t,UVRC.x,UVRC.u,1,[],UVRC.par);

%scrubber_gas V_001
V_001.Qg_in		= P_003.Qg;
V_001.T_in		= P_003.T;
V_001.Qg		= PV_001.Qg;
V_001.u = [...
    V_001.Qg_in; ...
    V_001.T_in; ...
    V_001.Qg];
V_001.dxdt = scrubber_gas(t,V_001.x,V_001.u,1,[],V_001.par);

%PI_pos LIC_001
LIC_001.dxdt = PI_pos(t,LIC_001.x,LIC_001.u,1,[],LIC_001.par);

%PI_pos LIC_002
LIC_002.dxdt = PI_pos(t,LIC_002.x,LIC_002.u,1,[],LIC_002.par);

%PI_pos LIC_003
LIC_003.dxdt = PI_pos(t,LIC_003.x,LIC_003.u,1,[],LIC_003.par);

%PI_pos LIC_004
LIC_004.dxdt = PI_pos(t,LIC_004.x,LIC_004.u,1,[],LIC_004.par);

%PI_pos LIC_005
LIC_005.dxdt = PI_pos(t,LIC_005.x,LIC_005.u,1,[],LIC_005.par);

%PI_pos LIC_006
LIC_006.dxdt = PI_pos(t,LIC_006.x,LIC_006.u,1,[],LIC_006.par);

%PI_pos PIC_001
PIC_001.dxdt = PI_pos(t,PIC_001.x,PIC_001.u,1,[],PIC_001.par);

%PI_pos TIC_001
TIC_001.dxdt = PI_pos(t,TIC_001.x,TIC_001.u,1,[],TIC_001.par);

%PI_pos TIC_002
TIC_002.dxdt = PI_pos(t,TIC_002.x,TIC_002.u,1,[],TIC_002.par);

%dxdt allocation:
dxdt=[...
    SG_001.dxdt; ...
	V_TO_001.dxdt; ...
	TO_001.dxdt; ...
	V_TO_002.dxdt; ...
	TO_002.dxdt; ...
	UVRC.dxdt; ...
	V_001.dxdt; ...
	LIC_001.dxdt; ...
	LIC_002.dxdt; ...
	LIC_003.dxdt; ...
	LIC_004.dxdt; ...
	LIC_005.dxdt; ...
	LIC_006.dxdt; ...
	PIC_001.dxdt; ...
	TIC_001.dxdt; ...
	TIC_002.dxdt; ...
	(PV_001.Qg)*1e-3/(24*3600)];

%----- State derivatives for sys output flag=1
if abs(imag(dxdt))>=eps
    error('Solution is imag');
end
sys=real(dxdt);
sys=sys(:);
%--------------------------------------------------------------------
    case 2,
sys=[];
%----- OUPUTS -------------------------------------------------------
    case 3,
%----- Output variables allocation:
y=[	SG_001.x; ...
	V_TO_001.x; ...
	TO_001.y(1); TO_001.y(2); TO_001.y(4); ...
	V_TO_002.x; ...
	TO_002.y(1); TO_002.y(2); TO_002.y(4); ...
	V_001.x; ...
	P_002.T; ...
	P_003.T; ...
	UVRC.scrb(1).T; UVRC.scrb(2).T; ...
	UVRC.comp(1).Qg; UVRC.comp(2).Qg; ...
	stk.Qf; stk.Qg; stk.BSW; stk.T; ...
	SG_001.Qw; ...
	TO_001.Qw; ...
	TO_002.Qw; ...
	SG_001.Qg; ...
	V_TO_001.Qg; ...
	V_TO_002.Qg];

%----- Outputs for sys output flag=3:
sys=y;
sys=real(sys(:));
%--------------------------------------------------------------------
%----- S-FUNCTION OTHER FLAGS ---------------------------------------
    case 4,
sys=[];
    case 9,
sys=[];
    otherwise   % Unexpected flags
DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
end
%--------------------------------------------------------------------    
end