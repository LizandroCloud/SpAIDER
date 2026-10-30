classdef oil_treater < ModelTypes  
    properties
        L 		= 12;      	%Vessel length          			[m]
        D 		= 4;    	%Vessel diameter        			[m]
        V 		= [];     	%Vessel total volume    			[m^3]
		TRL50 	= 300;		%Residence time for 50% BSW removal [s]
		b 		= 4;		%effWL 2nd parameter correlation 	[ ]
		TOGinf 	= 20;		%Minimum TOG achieved parameter 	[ml/m^3]
		TOGtau 	= 500;		%TOG residence time decay 			[s]
    end
    
    methods
        function obj = oil_treater(Name,varargin)
        %MODEL BRIEF
            obj.Description = 'Horizontal vessel for oil-water separation';
            obj.Name = Name;

        %PARAMETERS
        %User-defined parameters
            if length(varargin) > 12
                error('Too many input parameters as arguments. Valid parameters: {L,D,TRL50,b,TOGinf,TOGtau}');
            elseif length(varargin) > 1
                k=1;
                while k<length(varargin)
					switch varargin{k}
						case 'L'
							obj.L = varargin{k+1};
						case 'D'
							obj.D = varargin{k+1};
						case 'TRL50'
							obj.TRL50 = varargin{k+1};
						case 'b'
							obj.b = varargin{k+1};
						case 'TOGinf'
							obj.TOGinf = varargin{k+1};
						case 'TOGtau'
							obj.TOGtau = varargin{k+1};
					end
                    k=k+1;
                end
            end

        %Dependent parameters       
            obj.V = 0.25*pi*obj.L*obj.D^2;

        %Parameter tokens
            Lpar = num2str(obj.L);
            Dpar = num2str(obj.D);
            Vpar = num2str(obj.V);
			TRL50par = num2str(obj.TRL50);
			bpar = num2str(obj.b);
			TOGinfpar = num2str(obj.TOGinf);
			TOGtaupar = num2str(obj.TOGtau);
            
        %EQUIPMENT
            AddEquipment(obj,str3ph(strcat(obj.Name,'.Inlet')));
            AddEquipment(obj,liquid_stream(strcat(obj.Name,'.OutletL')));
            AddEquipment(obj,water_stream(strcat(obj.Name,'.OutletW')));

        %VARIABLES
            AddVariable(obj,'VOLstd',volume);
            AddVariable(obj,'VWLstd',volume);
            AddVariable(obj,'Mws',mass);
            AddVariable(obj,'VLstd',volume);
            AddVariable(obj,'VWstd',volume);
            AddVariable(obj,'E',energy);
            AddVariable(obj,'LevelL',distance);
			AddVariable(obj,'LevelW',distance);
			AddVariable(obj,'effWL',admensional);
			AddVariable(obj,'effOW',admensional);
			AddVariable(obj,'TRL',time_sec);
            AddVariable(obj,'dBodT',admensional);
            AddVariable(obj,'dBwdT',admensional);
            AddVariable(obj,'dEdT',admensional);

            idx = GetVarIndex(obj,'VOLstd');    obj.Variables{idx}.Name = strcat(Name,'.VOLstd');   obj.Variables{idx}.Description = 'Oil std. vol. hold-up at liquid phase';   	obj.Variables{idx}.Value = 150;		obj.Variables{idx}.Unit = 'Sm^3';
            idx = GetVarIndex(obj,'VWLstd');    obj.Variables{idx}.Name = strcat(Name,'.VWLstd');   obj.Variables{idx}.Description = 'Water std. vol. hold-up at liquid phase'; 	obj.Variables{idx}.Value = 1.00;  	obj.Variables{idx}.Unit = 'Sm^3';
            idx = GetVarIndex(obj,'Mws');       obj.Variables{idx}.Name = strcat(Name,'.Mws');      obj.Variables{idx}.Description = 'Total water salt mass hold-up';           	obj.Variables{idx}.Value = 1e5; 	obj.Variables{idx}.Unit = 'kg';
            idx = GetVarIndex(obj,'VLstd');     obj.Variables{idx}.Name = strcat(Name,'.VLstd');    obj.Variables{idx}.Description = 'Liquid phase std. vol. hold-up';          	obj.Variables{idx}.Value = 150; 	obj.Variables{idx}.Unit = 'Sm^3';
            idx = GetVarIndex(obj,'VWstd');     obj.Variables{idx}.Name = strcat(Name,'.VWstd');    obj.Variables{idx}.Description = 'Water phase std. vol. hold-up';           	obj.Variables{idx}.Value = 150; 	obj.Variables{idx}.Unit = 'Sm^3';
            idx = GetVarIndex(obj,'E');         obj.Variables{idx}.Name = strcat(Name,'.E');        obj.Variables{idx}.Description = 'Energy hold-up';                          	obj.Variables{idx}.Value = 5e6; 	obj.Variables{idx}.Unit = 'kJ';
            idx = GetVarIndex(obj,'LevelL');    obj.Variables{idx}.Name = strcat(Name,'.LevelL');   obj.Variables{idx}.Description = 'Liquid level';                            	obj.Variables{idx}.Value = 1.2; 	obj.Variables{idx}.Unit = 'm';
			idx = GetVarIndex(obj,'LevelW');    obj.Variables{idx}.Name = strcat(Name,'.LevelW');   obj.Variables{idx}.Description = 'Water level';                            		obj.Variables{idx}.Value = 1.2; 	obj.Variables{idx}.Unit = 'm';
			idx = GetVarIndex(obj,'effWL');    	obj.Variables{idx}.Name = strcat(Name,'.effWL');   	obj.Variables{idx}.Description = 'Water removal efficiency from light phase';   obj.Variables{idx}.Value = 0.2; 	obj.Variables{idx}.Unit = '';
			idx = GetVarIndex(obj,'effOW');    	obj.Variables{idx}.Name = strcat(Name,'.effOW');   	obj.Variables{idx}.Description = 'Oil removal efficiency from water phase';   	obj.Variables{idx}.Value = 0.999;   obj.Variables{idx}.Unit = '';
			idx = GetVarIndex(obj,'TRL');    	obj.Variables{idx}.Name = strcat(Name,'.TRL');   	obj.Variables{idx}.Description = 'Light phase residence time';   				obj.Variables{idx}.Value = 300; 	obj.Variables{idx}.Unit = 's';
            idx = GetVarIndex(obj,'dBodT');    	obj.Variables{idx}.Name = strcat(Name,'.dBodT');   	obj.Variables{idx}.Description = 'Oil FVF derivative to T';                     obj.Variables{idx}.Value = 1.5e-3; 	obj.Variables{idx}.Unit = '1/K';
            idx = GetVarIndex(obj,'dBwdT');    	obj.Variables{idx}.Name = strcat(Name,'.dBwdT');   	obj.Variables{idx}.Description = 'Water FVF derivative to T';                   obj.Variables{idx}.Value = -6e-4; 	obj.Variables{idx}.Unit = '1/K';
            idx = GetVarIndex(obj,'dEdT');    	obj.Variables{idx}.Name = strcat(Name,'.dEdT');   	obj.Variables{idx}.Description = 'Energy hold-up derivative to T';              obj.Variables{idx}.Value = 5e5; 	obj.Variables{idx}.Unit = 'kJ/K';
			

        %Variable tokens
            VOLstd              = strcat(obj.Name,'.VOLstd');
            VWLstd              = strcat(obj.Name,'.VWLstd');
            Mws                 = strcat(obj.Name,'.Mws');
            VLstd               = strcat(obj.Name,'.VLstd');
            VWstd               = strcat(obj.Name,'.VWstd');
			E                   = strcat(obj.Name,'.E');
            LevelL              = strcat(obj.Name,'.LevelL');
			LevelW              = strcat(obj.Name,'.LevelW');
			effWL               = strcat(obj.Name,'.effWL');
			effOW               = strcat(obj.Name,'.effOW');
			TRL               	= strcat(obj.Name,'.TRL');
            dBodT               = strcat(obj.Name,'.dBodT');
            dBwdT               = strcat(obj.Name,'.dBwdT');
            dEdT               	= strcat(obj.Name,'.dEdT');
            
            Inlet.Qf            = strcat(obj.Name,'.Inlet.Qf');
            Inlet.Qg            = strcat(obj.Name,'.Inlet.Qg');
            Inlet.BSW           = strcat(obj.Name,'.Inlet.BSW');
            Inlet.P             = strcat(obj.Name,'.Inlet.P');
            Inlet.T             = strcat(obj.Name,'.Inlet.T');
            Inlet.eO            = strcat(obj.Name,'.Inlet.eO');
            Inlet.eG            = strcat(obj.Name,'.Inlet.eG');
            Inlet.eW            = strcat(obj.Name,'.Inlet.eW');
            Inlet.rho_Wstd      = strcat(obj.Name,'.Inlet.rho_Wstd');
            Inlet.Rs            = strcat(obj.Name,'.Inlet.Rs');
            Inlet.Rv            = strcat(obj.Name,'.Inlet.Rv');
            Inlet.Bo            = strcat(obj.Name,'.Inlet.Bo');
            Inlet.Bw            = strcat(obj.Name,'.Inlet.Bw');
            Inlet.d60           = strcat(obj.Name,'.Inlet.d60');
            Inlet.d60G          = strcat(obj.Name,'.Inlet.d60G');
            Inlet.Tbn           = strcat(obj.Name,'.Inlet.Tbn');
            Inlet.TbnG          = strcat(obj.Name,'.Inlet.TbnG');
            Inlet.dg            = strcat(obj.Name,'.Inlet.dg');
            Inlet.dgs           = strcat(obj.Name,'.Inlet.dgs');
            Inlet.Ycg           = strcat(obj.Name,'.Inlet.Ycg');
            Inlet.Ycgs          = strcat(obj.Name,'.Inlet.Ycgs');
            Inlet.ws            = strcat(obj.Name,'.Inlet.ws');

            OutletL.Qf          = strcat(obj.Name,'.OutletL.Qf');
            OutletL.P           = strcat(obj.Name,'.OutletL.P');
            OutletL.T           = strcat(obj.Name,'.OutletL.T');
            OutletL.BSW         = strcat(obj.Name,'.OutletL.BSW');
            OutletL.eO          = strcat(obj.Name,'.OutletL.eO');
            OutletL.eW          = strcat(obj.Name,'.OutletL.eW');
            OutletL.rho_Wstd    = strcat(obj.Name,'.OutletL.rho_Wstd');
            OutletL.Rs          = strcat(obj.Name,'.OutletL.Rs');
            OutletL.Bo          = strcat(obj.Name,'.OutletL.Bo');
            OutletL.Bw          = strcat(obj.Name,'.OutletL.Bw');
            OutletL.d60         = strcat(obj.Name,'.OutletL.d60');
            OutletL.Tbn         = strcat(obj.Name,'.OutletL.Tbn');
            OutletL.dg          = strcat(obj.Name,'.OutletL.dg');
            OutletL.Ycg         = strcat(obj.Name,'.OutletL.Ycg');
            OutletL.ws          = strcat(obj.Name,'.OutletL.ws');

            OutletW.Qw          = strcat(obj.Name,'.OutletW.Qw');
            OutletW.P           = strcat(obj.Name,'.OutletW.P');
            OutletW.T           = strcat(obj.Name,'.OutletW.T');
            OutletW.TOG         = strcat(obj.Name,'.OutletW.TOG');
            OutletW.eW          = strcat(obj.Name,'.OutletW.eW');
			OutletW.rhoW    	= strcat(obj.Name,'.OutletW.rhoW');
            OutletW.rho_Wstd    = strcat(obj.Name,'.OutletW.rho_Wstd');
            OutletW.d60         = strcat(obj.Name,'.OutletW.d60');
            OutletW.Tbn         = strcat(obj.Name,'.OutletW.Tbn');
            OutletW.ws          = strcat(obj.Name,'.OutletW.ws');          
            
        %REVERSE MODEL INHERITANCE
            SystemVariables(obj);
            SystemEquations(obj);
            SystemInitConditions(obj);
        
        %EQUATIONS
        %Gas std. vol. balance
            AddEquation(obj,strcat(OutletL.Rs,' - ',Inlet.Rs));
            
        %Water std. vol. balance
            AddEquation(obj,strcat(...
                'diff(',VWstd,')',' - (',Inlet.Qf,'/(24*3600)*',Inlet.BSW,'*',effWL,' - ',OutletW.Qw,'/(24*3600)',')'));

        %Water salt mass balance
            AddEquation(obj,strcat(...
                'diff(',Mws,')',' - (', ...
                Inlet.Qf,'/(24*3600)*',Inlet.BSW,'*',effWL,'*',Inlet.rho_Wstd,'*',Inlet.ws, ...
                    ' - ',OutletW.Qw,'/(24*3600)*',OutletW.rho_Wstd,'*',OutletW.ws,...
                ')'));

        %Oil std. vol. hold-up at liquid phase
            AddEquation(obj,strcat(VOLstd,' - ',VLstd,'*(1-',OutletL.BSW,')'));

        %Water std. vol. hold-up at liquid phase
            AddEquation(obj,strcat(VWLstd,' - ',VLstd,'*',OutletL.BSW));
        
        %Water salt mass hold-up
            AddEquation(obj,strcat(Mws,' - ',VWstd,'*',OutletW.rho_Wstd,'*',OutletW.ws));

        %Water salt content constraint
            AddEquation(obj,strcat(OutletL.ws,' - ',OutletW.ws));

		%Height constraint
            AddEquation(obj,strcat(Dpar,' - (',LevelW,' + ',LevelL,')'));

        %Volume constraint
            AddEquation(obj,strcat(...
                Vpar,' - (',VWstd,'*',OutletL.Bw,' + ',VOLstd,'*',OutletL.Bo,' + ',VWLstd,'*',OutletL.Bw,')'));

%            AddEquation(obj,strcat(...
%                Vpar,' - (',VWstd,'*',OutletL.Bw,' + ',VLstd,'*((1-',OutletL.BSW,')*',OutletL.Bo,' + ',OutletL.BSW,'*',OutletL.Bw,'))'));

		%Water phase volume hold-up
            AddEquation(obj,strcat(...
                VWstd,'*',OutletL.Bw,' - (',...
                '(0.25*',Lpar,'*',Dpar,'^2)*(','acos(1-2*',LevelW,'/',Dpar,')',...
                    '- (1-2*',LevelW,'/',Dpar,')*sin(acos(1-2*',LevelW,'/',Dpar,')))',...
                ')'));

		%Relation between BSW and residence time for light phase
            AddEquation(obj,strcat(...
                OutletL.BSW,' - ',Inlet.BSW,'/(1+(',TRL,'/',TRL50par,')^',bpar,')'));

%            AddEquation(obj,strcat(...
%                OutletL.BSW,'*','(1+(',TRL,'/',TRL50par,')^',bpar,')',' - ',Inlet.BSW));

		%Oil residence time at light phase
            AddEquation(obj,strcat(...
                VOLstd,'*',OutletL.Bo,' - ',TRL,'*',Inlet.Qf,'/(24*3600)*(1-',Inlet.BSW,')*',Inlet.Bo));

%            AddEquation(obj,strcat(...
%                VLstd,'*(1-',OutletL.BSW,')*',OutletL.Bo,' - ',TRL,'*',Inlet.Qf,'/(24*3600)*(1-',Inlet.BSW,')*',Inlet.Bo));

		%Water removal efficiency from light phase
            AddEquation(obj,strcat(...
                effWL,'*',Inlet.Qf,'/(24*3600)*',Inlet.BSW,' - (', ...
					Inlet.Qf,'/(24*3600)*',Inlet.BSW,' - ',OutletL.Qf,'/(24*3600)*',OutletL.BSW,')'));

		%TOG correlation for water outlet stream
            AddEquation(obj,strcat(OutletW.TOG,' - (',TOGinfpar,' + ','(1000 - ',TOGinfpar,')*', ...
					'exp(-((',VWstd,'*',OutletL.Bw,')/(',Inlet.Qf,'/(24*3600)*',Inlet.BSW,'*',Inlet.Bw,'))/',TOGtaupar,')', ...
				')'));

%            AddEquation(obj,strcat('(',OutletW.TOG,' - ',TOGinfpar,')/(','1000 - ',TOGinfpar,') - (', ...
%					'exp(-((',VWstd,'*',OutletL.Bw,')/(',Inlet.Qf,'/(24*3600)*',Inlet.BSW,'*',Inlet.Bw,'))/',TOGtaupar,')', ...
%				')'));

%            AddEquation(obj,strcat(OutletW.TOG,' - 0.00'));

		%Oil removal efficiency from water phase
            AddEquation(obj,strcat(...
                effOW,'*',Inlet.Qf,'/(24*3600)*(1-',Inlet.BSW,')',' - (', ...
					Inlet.Qf,'/(24*3600)*(1-',Inlet.BSW,')',' - ',OutletW.Qw,'/(24*3600)*',OutletW.TOG,'*1e-6',')'));

%            AddEquation(obj,strcat(effOW,' - ', ...
%					'(',Inlet.Qf,'/(24*3600)*(1-',Inlet.BSW,')',' - ',OutletW.Qw,'/(24*3600)*',OutletW.TOG,'*1e-6',')/(',Inlet.Qf,'/(24*3600)*(1-',Inlet.BSW,')',')'));

%             AddEquation(obj,strcat(effOW,' - 1.00'));

        %Energy balance
            AddEquation(obj,strcat(...
                'diff(',E,')',' - (', ...
                Inlet.Qf,'/(24*3600)*((1-',Inlet.BSW,')*',effOW,'*',Inlet.d60,'*PP.rho_H2Ostd*',Inlet.eO, ...
                    ' + ',Inlet.BSW,'*',Inlet.rho_Wstd,'*',Inlet.eW,')', ...
                ' - ',OutletL.Qf,'/(24*3600)*((1-',OutletL.BSW,')*',OutletL.d60,'*PP.rho_H2Ostd*',OutletL.eO, ...
                    ' + ',OutletL.BSW,'*',OutletL.rho_Wstd,'*',OutletL.eW,')', ...
                ' - ',OutletW.Qw,'/(24*3600)*',OutletW.rho_Wstd,'*',OutletW.eW,...
                ')'));

%            AddEquation(obj,strcat(OutletL.T,' - ',Inlet.T));

        %Energy hold-up
            AddEquation(obj,strcat(...
                E,' - (',VOLstd,'*',OutletL.d60,'*PP.rho_H2Ostd*',OutletL.eO, ...
                    ' + ',VWLstd,'*',OutletL.rho_Wstd,'*',OutletL.eW,' + ',VWstd,'*',OutletW.rho_Wstd,'*',OutletW.eW,' - ',OutletL.P,'*',Vpar,...
                ')'));

%            AddEquation(obj,strcat(...
%                E,' - (',VLstd,'*((1-',OutletL.BSW,')*',OutletL.d60,'*PP.rho_H2Ostd*',OutletL.eO,' + ',OutletL.BSW,'*',OutletL.rho_Wstd,'*',OutletL.eW,')', ...
%                    ' + ',VWstd,'*',OutletW.rho_Wstd,'*',OutletW.eW,' - ',OutletL.P,'*',Vpar,...
%                ')'));

        %Volume change rate constraint (degasser liquid outlet gravity flow constraint)
%            AddEquation(obj,strcat(...
%				Inlet.Qf,'/(24*3600)*(',OutletL.Bo,'*(1-',Inlet.BSW,')*',effOW,' + ',OutletL.Bw,'*',Inlet.BSW,') + diff(',E,')*', ...
%                    '(',VOLstd,'*',dBodT,'/',dEdT,' + ',VWLstd,'*',dBwdT,'/',dEdT,' + ',VWstd,'*',dBwdT,'/',dEdT,')',' - (', ... 
%				OutletL.Qf,'/(24*3600)*(',OutletL.Bo,'*(1-',OutletL.BSW,') + ',OutletL.Bw,'*',OutletL.BSW,') + ',OutletW.Qw,'/(24*3600)*',OutletL.Bw, ...
%                ')'));

            AddEquation(obj,strcat(...
				Inlet.Qf,'/(24*3600)*(',OutletL.Bo,'*(1-',Inlet.BSW,')*',effOW,' + ',OutletL.Bw,'*',Inlet.BSW,') + diff(',E,')/',dEdT,'*', ...
                    '(',VOLstd,'*',dBodT,' + ',VWLstd,'*',dBwdT,' + ',VWstd,'*',dBwdT,')',' - (', ... 
				OutletL.Qf,'/(24*3600)*(',OutletL.Bo,'*(1-',OutletL.BSW,') + ',OutletL.Bw,'*',OutletL.BSW,') + ',OutletW.Qw,'/(24*3600)*',OutletL.Bw, ...
                ')'));

%            AddEquation(obj,strcat(...
%				Inlet.Qf,'/(24*3600)*(',OutletL.Bo,'*(1-',Inlet.BSW,')*',effOW,' + ',OutletL.Bw,'*',Inlet.BSW,') + diff(',E,')/',dEdT,'*', ...
%                    '(',VLstd,'*((1-',OutletL.BSW,')*',dBodT,' + ',OutletL.BSW,'*',dBwdT,') + ',VWstd,'*',dBwdT,')',' - (', ... 
%				OutletL.Qf,'/(24*3600)*(',OutletL.Bo,'*(1-',OutletL.BSW,') + ',OutletL.Bw,'*',OutletL.BSW,') + ',OutletW.Qw,'/(24*3600)*',OutletL.Bw, ...
%                ')'));

%            AddEquation(obj,strcat(...
%				Inlet.Qf,'/(24*3600)*(',OutletL.Bo,'*(1-',Inlet.BSW,')*',effOW,' + ',OutletL.Bw,'*',Inlet.BSW,')',' - (', ... 
%				OutletL.Qf,'/(24*3600)*(',OutletL.Bo,'*(1-',OutletL.BSW,') + ',OutletL.Bw,'*',OutletL.BSW,') + ',OutletW.Qw,'/(24*3600)*',OutletL.Bw, ...
%                ')'));

%            AddEquation(obj,strcat(Inlet.Qf,'/(24*3600)*(1-',Inlet.BSW,')*',effOW,' - ',OutletL.Qf,'/(24*3600)*(1-',OutletL.BSW,')'));

		%dBo/dT
			AddEquation(obj,strcat(dBodT,' - dOilFVFdT(PP,',OutletL.Rs,',',OutletL.d60,',',OutletL.dg,')'));
%            AddEquation(obj,strcat(dBodT,' - 9.5e-4'));

		%dBw/dT
			AddEquation(obj,strcat(dBwdT,' - dWaterFVFdT(PP,',OutletL.P,',',OutletL.T,')'));
%            AddEquation(obj,strcat(dBwdT,' - -6.e-4'));

		%dE/dT
            AddEquation(obj,strcat(...
                dEdT,' - (',VOLstd,'*',OutletL.d60,'*PP.rho_H2Ostd*OilSpecHeatCap(PP,',OutletL.d60,',',OutletL.Tbn,')', ...
                    ' + (',VWLstd,' + ',VWstd,')*',OutletW.rho_Wstd,'*WaterSpecHeatCap(PP,',OutletL.P,',',OutletL.T,',',OutletL.ws,')',...
                ')'));

%            AddEquation(obj,strcat(...
%                dEdT,' - (',VLstd,'*(1-',OutletL.BSW,')*',OutletL.d60,'*PP.rho_H2Ostd*OilSpecHeatCap(PP,',OutletL.d60,',',OutletL.Tbn,')', ...
%                    ' + (',VLstd,'*',OutletL.BSW,' + ',VWstd,')*',OutletW.rho_Wstd,'*WaterSpecHeatCap(PP,',OutletL.P,',',OutletL.T,',',OutletL.ws,')',...
%                ')'));

%            AddEquation(obj,strcat(dEdT,' - 5e5'));


%            AddEquation(obj,strcat(...
%                dEdT,' - (',VOLstd,'*',OutletL.d60,'*PP.rho_H2Ostd*OilSpecHeatCap(PP,',OutletL.d60,',',OutletL.Tbn,')', ...
%                    ' + (',VWLstd,' + ',VWstd,')*',OutletW.rho_Wstd,'*(4.180 + ',OutletW.rhoW,'*',OutletW.ws,'*(-4.396e-3 + 4.8e-6*',OutletW.rhoW,'*',OutletW.ws,'))',...
%                ')'));

        %Property mixing rule for oil specific gravity (std. vol. balance)
            AddEquation(obj,strcat(...
                OutletL.d60,'*diff(',VOLstd,') + ',VOLstd,'*diff(',OutletL.d60,')',' - (', ...
                Inlet.Qf,'/(24*3600)*(1-',Inlet.BSW,')*',effOW,'*',Inlet.d60,...
                    ' - ',OutletL.Qf,'/(24*3600)*(1-',OutletL.BSW,')*',OutletL.d60, ...
                ')'));

%            AddEquation(obj,strcat(...
%                '(1-',OutletL.BSW,')*',OutletL.d60,'*diff(',VLstd,') - ',VLstd,'*',OutletL.d60,'*diff(',OutletL.BSW,') + ',VLstd,'*(1-',OutletL.BSW,')*diff(',OutletL.d60,')',' - (', ...
%                Inlet.Qf,'/(24*3600)*(1-',Inlet.BSW,')*',effOW,'*',Inlet.d60,...
%                    ' - ',OutletL.Qf,'/(24*3600)*(1-',OutletL.BSW,')*',OutletL.d60, ...
%                ')'));

%            AddEquation(obj,strcat(OutletL.d60,' - ',Inlet.d60));

        %Property mixing rule for Tbn by means of Watson K-factor (std. vol. balance)
            AddEquation(obj,strcat(...
                '(',OutletL.Tbn,')^(1/3)','*diff(',VOLstd,') + ',VOLstd,'*1/3*',OutletL.Tbn,'^(-2/3)*diff(',OutletL.Tbn,')',' - (', ...
                	Inlet.Qf,'/(24*3600)*(1-',Inlet.BSW,')*',effOW,'*(',Inlet.Tbn,')^(1/3)', ...
                    ' - ',OutletL.Qf,'/(24*3600)*(1-',OutletL.BSW,')*(',OutletL.Tbn,')^(1/3)',...
                ')'));

%            AddEquation(obj,strcat(...
%                '(1-',OutletL.BSW,')*(',OutletL.Tbn,')^(1/3)*diff(',VLstd,') - ',VLstd,'*(',OutletL.Tbn,')^(1/3)*diff(',OutletL.BSW,') + ', ...
%                    VLstd,'*(1-',OutletL.BSW,')*(1/3)*(',OutletL.Tbn,')^(-2/3)*diff(',OutletL.Tbn,')',' - (', ...
%                Inlet.Qf,'/(24*3600)*(1-',Inlet.BSW,')*',effOW,'*(',Inlet.Tbn,')^(1/3)', ...
%                    ' - ',OutletL.Qf,'/(24*3600)*(1-',OutletL.BSW,')*(',OutletL.Tbn,')^(1/3)',...
%                ')'));
            
%            AddEquation(obj,strcat(OutletL.Tbn,' - ',Inlet.Tbn));

        %Equates of liquid outlet stream compositional variables
            AddEquation(obj,strcat(OutletL.dg,' - ',Inlet.dg));
			AddEquation(obj,strcat(OutletL.Ycg,' - ',Inlet.Ycg));

        %Equates of water outlet stream compositional variables
            AddEquation(obj,strcat(OutletW.d60,' - ',OutletL.d60));
			AddEquation(obj,strcat(OutletW.Tbn,' - ',OutletL.Tbn));

        %Thermal equlibrium
            AddEquation(obj,strcat(OutletL.T,' - ',OutletW.T));

        %Mechanical equlibrium
            AddEquation(obj,strcat(Inlet.P,' - ',OutletL.P));
            AddEquation(obj,strcat(OutletW.P,' - ',OutletL.P));

        end
    end
    
end

