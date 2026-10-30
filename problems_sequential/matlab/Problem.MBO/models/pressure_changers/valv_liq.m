classdef valv_liq < ModelTypes  
    properties
        Cv 		= 100;      %Valve coefficient constant US units	[gal/min]
        Kv 		= [];    	%Valve coefficient constant SI units    [m^3/h]
        Av 		= [];     	%Valve coefficient constant    			[m^2]
		FL 		= 0.91;		%Liquid pressure recovery factor 		[]	
	end
	
	methods
        function obj = valv_liq(Name,varargin)
        %MODEL BRIEF
            obj.Description = 'Basic equations for liquid stream flowing through valves';
            obj.Name = Name;

        %PARAMETERS
        %User-defined parameters
            if length(varargin) > 6
                error('Too many input parameters as arguments. Valid parameters: {Cv,FL,tau}');
            elseif length(varargin) > 1
                k=1;
                while k<length(varargin)
					switch varargin{k}
						case 'Cv'
							obj.Cv = varargin{k+1};
						case 'FL'
							obj.FL = varargin{k+1};
					end
                    k=k+1;
                end
            end

        %Dependent parameters       
            obj.Kv = obj.Cv*(0.227124707)*sqrt(14.696/1.01325);
			obj.Av = obj.Cv*(0.0000630902)*sqrt(999.054*14.696/101325);

        %Parameter tokens
            Avpar = num2str(obj.Av);
			FLpar = num2str(obj.FL);
            
        %EQUIPMENT
            AddEquipment(obj,str3ph(strcat(obj.Name,'.Inlet')));
            AddEquipment(obj,triphase_stream(strcat(obj.Name,'.Outlet')));

        %VARIABLES
            AddVariable(obj,'fc',admensional);
			AddVariable(obj,'c',admensional);
			AddVariable(obj,'FF',admensional);
			AddVariable(obj,'DP',pressure);
			AddVariable(obj,'DPchok',pressure);
			AddVariable(obj,'PvL',pressure);
			AddVariable(obj,'PpcL',pressure);
			AddVariable(obj,'PpcO',pressure);
			AddVariable(obj,'MwO',molweight);
            AddVariable(obj,'Rseq',admensional);
            AddVariable(obj,'SumInletGas',flow_vol);
            AddVariable(obj,'Gunbal',flow_vol);

			idx = GetVarIndex(obj,'fc');    		obj.Variables{idx}.Name = strcat(Name,'.fc');   		obj.Variables{idx}.Description = 'Valve inherent characteristic function';  obj.Variables{idx}.Value = 0.5; 	obj.Variables{idx}.Unit = '';
			idx = GetVarIndex(obj,'c');    			obj.Variables{idx}.Name = strcat(Name,'.c');	   		obj.Variables{idx}.Description = 'Valve stem position';  					obj.Variables{idx}.Value = 0.5; 	obj.Variables{idx}.Unit = '';
			idx = GetVarIndex(obj,'FF');    		obj.Variables{idx}.Name = strcat(Name,'.FF');	   		obj.Variables{idx}.Description = 'Liquid critical pressure ratio factor';  	obj.Variables{idx}.Value = 0.96;	obj.Variables{idx}.Unit = '';
            idx = GetVarIndex(obj,'DP');    		obj.Variables{idx}.Name = strcat(Name,'.DP');   		obj.Variables{idx}.Description = 'Pressure drop';   						obj.Variables{idx}.Value = 10;  	obj.Variables{idx}.Unit = 'kPa';
			idx = GetVarIndex(obj,'DPchok');    	obj.Variables{idx}.Name = strcat(Name,'.DPchok');   	obj.Variables{idx}.Description = 'Pressure drop at choked flow';   			obj.Variables{idx}.Value = 5;  		obj.Variables{idx}.Unit = 'kPa';
            idx = GetVarIndex(obj,'PvL');    		obj.Variables{idx}.Name = strcat(Name,'.PvL');   		obj.Variables{idx}.Description = 'Liquid vapour pressure';   				obj.Variables{idx}.Value = 80;  	obj.Variables{idx}.Unit = 'kPa';
			idx = GetVarIndex(obj,'PpcL');    		obj.Variables{idx}.Name = strcat(Name,'.PpcL');   		obj.Variables{idx}.Description = 'Liquid pseudo-critical pressure';   		obj.Variables{idx}.Value = 30;  	obj.Variables{idx}.Unit = 'kPa';
			idx = GetVarIndex(obj,'PpcO');    		obj.Variables{idx}.Name = strcat(Name,'.PpcO');   		obj.Variables{idx}.Description = 'Oil pseudo-critical pressure';   			obj.Variables{idx}.Value = 300;  	obj.Variables{idx}.Unit = 'kPa';
			idx = GetVarIndex(obj,'MwO');    		obj.Variables{idx}.Name = strcat(Name,'.MwO');   		obj.Variables{idx}.Description = 'Molar mass for oil phase';   				obj.Variables{idx}.Value = 300;  	obj.Variables{idx}.Unit = 'kg/kmol';
			idx = GetVarIndex(obj,'Rseq');    		obj.Variables{idx}.Name = strcat(Name,'.Rseq');   		obj.Variables{idx}.Description = 'Solution gas-oil ratio at eq.';      		obj.Variables{idx}.Value = 10;  	obj.Variables{idx}.Unit = 'Sm^3/Sm^3';
            idx = GetVarIndex(obj,'SumInletGas');   obj.Variables{idx}.Name = strcat(Name,'.SumInletGas');  obj.Variables{idx}.Description = 'Inlet streams gas std. vol. flow sum';   	obj.Variables{idx}.Value = 2e5;		obj.Variables{idx}.Unit = 'Sm^3/d';
            idx = GetVarIndex(obj,'Gunbal');    	obj.Variables{idx}.Name = strcat(Name,'.Gunbal');   	obj.Variables{idx}.Description = 'Gas std. vol. unbalance'; 				obj.Variables{idx}.Value = 2e5; 	obj.Variables{idx}.Unit = 'Sm^3/d';

        %Variable tokens
            fc          		= strcat(obj.Name,'.fc');
			c          			= strcat(obj.Name,'.c');
			FF          		= strcat(obj.Name,'.FF');
			DP          		= strcat(obj.Name,'.DP');
            DPchok      		= strcat(obj.Name,'.DPchok');
			PvL      			= strcat(obj.Name,'.PvL');
			PpcL      			= strcat(obj.Name,'.PpcL');
			PpcO      			= strcat(obj.Name,'.PpcO');
			MwO      			= strcat(obj.Name,'.MwO');
			Rseq        		= strcat(obj.Name,'.Rseq');
            SumInletGas			= strcat(obj.Name,'.SumInletGas');
            Gunbal      		= strcat(obj.Name,'.Gunbal');

            Inlet.Qf            = strcat(obj.Name,'.Inlet.Qf');
            Inlet.Qg            = strcat(obj.Name,'.Inlet.Qg'); 
            Inlet.BSW           = strcat(obj.Name,'.Inlet.BSW');
            Inlet.P             = strcat(obj.Name,'.Inlet.P');
            Inlet.T             = strcat(obj.Name,'.Inlet.T');
            Inlet.eO            = strcat(obj.Name,'.Inlet.eO');
            Inlet.eG            = strcat(obj.Name,'.Inlet.eG');
            Inlet.eW            = strcat(obj.Name,'.Inlet.eW');
            Inlet.rhoL          = strcat(obj.Name,'.Inlet.rhoL');
            Inlet.rho_Wstd      = strcat(obj.Name,'.Inlet.rho_Wstd');
            Inlet.Rs            = strcat(obj.Name,'.Inlet.Rs');
            Inlet.Rv            = strcat(obj.Name,'.Inlet.Rv');
            Inlet.d60           = strcat(obj.Name,'.Inlet.d60');
            Inlet.d60G          = strcat(obj.Name,'.Inlet.d60G');
            Inlet.Tbn           = strcat(obj.Name,'.Inlet.Tbn');
            Inlet.TbnG          = strcat(obj.Name,'.Inlet.TbnG');
            Inlet.dg            = strcat(obj.Name,'.Inlet.dg');
            Inlet.dgs           = strcat(obj.Name,'.Inlet.dgs');
            Inlet.Ycg           = strcat(obj.Name,'.Inlet.Ycg');
            Inlet.Ycgs          = strcat(obj.Name,'.Inlet.Ycgs');
            Inlet.ws            = strcat(obj.Name,'.Inlet.ws');

            Outlet.Qf           = strcat(obj.Name,'.Outlet.Qf');
            Outlet.Qg           = strcat(obj.Name,'.Outlet.Qg'); 
            Outlet.BSW          = strcat(obj.Name,'.Outlet.BSW');
            Outlet.P            = strcat(obj.Name,'.Outlet.P');
            Outlet.T            = strcat(obj.Name,'.Outlet.T');
            Outlet.eO           = strcat(obj.Name,'.Outlet.eO');
            Outlet.eG           = strcat(obj.Name,'.Outlet.eG');
            Outlet.eW           = strcat(obj.Name,'.Outlet.eW');
            Outlet.rhoL         = strcat(obj.Name,'.Outlet.rhoL');
            Outlet.rho_Wstd     = strcat(obj.Name,'.Outlet.rho_Wstd');
            Outlet.Rs           = strcat(obj.Name,'.Outlet.Rs');
            Outlet.Rv           = strcat(obj.Name,'.Outlet.Rv');
            Outlet.d60          = strcat(obj.Name,'.Outlet.d60');
            Outlet.d60G         = strcat(obj.Name,'.Outlet.d60G');
            Outlet.Tbn          = strcat(obj.Name,'.Outlet.Tbn');
            Outlet.TbnG         = strcat(obj.Name,'.Outlet.TbnG');
            Outlet.dg           = strcat(obj.Name,'.Outlet.dg');
            Outlet.dgs          = strcat(obj.Name,'.Outlet.dgs');
            Outlet.Ycg          = strcat(obj.Name,'.Outlet.Ycg');
            Outlet.Ycgs         = strcat(obj.Name,'.Outlet.Ycgs');
            Outlet.ws           = strcat(obj.Name,'.Outlet.ws');     
            
        %REVERSE MODEL INHERITANCE
            SystemVariables(obj);
            SystemEquations(obj);
            SystemInitConditions(obj);
        
        %EQUATIONS
		%Pressure-flow ralation
            AddEquation(obj,strcat(...
                Inlet.Qf,'/(24*3600)*((1-',Inlet.BSW,')*',Inlet.d60,'*PP.rho_H2Ostd + ',Inlet.BSW,'*',Inlet.rho_Wstd,') - (', ...
					'0.5*(1+tanh(((',DPchok,'-',DP,')+1e-3)/1e-4))*',fc,'*',Avpar,'*sqrt(',Inlet.rhoL,'*',DP,'*1e3)', ...
					' + 0.5*(1+tanh((',DP,'-',DPchok,')/1e-4))*',fc,'*',Avpar,'*sqrt(',Inlet.rhoL,'*',DPchok,'*1e3)', ...
               ')'));
%            AddEquation(obj,strcat(...
%                '(',Inlet.Qf,'/(24*3600)*((1-',Inlet.BSW,')*',Inlet.d60,'*PP.rho_H2Ostd + ',Inlet.BSW,'*',Inlet.rho_Wstd,')/', ...
%                    '(',fc,'*',Avpar,'))^2 - (', ...
%					'0.5*(1+tanh(((',DPchok,'-',DP,')+1e-3)/1e-4))*',Inlet.rhoL,'*',DP,'*1e3', ...
%					'+ 0.5*(1+tanh((',DP,'-',DPchok,')/1e-4))*',Inlet.rhoL,'*',DPchok,'*1e3', ...
%               ')'));

        %Oil std. vol. balance
            AddEquation(obj,strcat(...
                Inlet.Qf,'/(24*3600)*(1-',Inlet.BSW,')',' - (', ...
					Outlet.Qf,'/(24*3600)*(1-',Outlet.BSW,') + ',Outlet.Qg,'/(24*3600)*',Outlet.Rv,'*1e-3', ...
                ')'));

        %Sum of inlet streams gas std. vol.
            AddEquation(obj,strcat( ...
				SumInletGas,'/(24*3600)',' - ',Inlet.Qf,'/(24*3600)*(1-',Inlet.BSW,')*',Inlet.Rs));

        %Gas std. vol. unbalance
            AddEquation(obj,strcat(...
                Gunbal,'/(24*3600)',' - (',Outlet.Qf,'/(24*3600)*(1-',Outlet.BSW,')*',Rseq,' - ',SumInletGas,'/(24*3600)',')'));

		%Std. vol. balance for liberated gas
            AddEquation(obj,strcat(...
                Outlet.Qg,'/(24*3600)',' - 0.5*(1+tanh(-',Gunbal,'/1e-6))*(-',Gunbal,'/(24*3600)',')'));

        %Std. vol. balance for solution gas-oil ratio
            AddEquation(obj,strcat(...
                Outlet.Qf,'/(24*3600)*(1-',Outlet.BSW,')*',Outlet.Rs,' - (', ...
					'0.5*(1+tanh(',Gunbal,'/1e-6))*',SumInletGas,'/(24*3600)', ...
					' + 0.5*(1+tanh(-',Gunbal,'/1e-6))*',Outlet.Qf,'/(24*3600)*(1-',Outlet.BSW,')*',Rseq, ...
				')'));
            
        %Water std. vol. balance
            AddEquation(obj,strcat(Inlet.Qf,'/(24*3600)*',Inlet.BSW,' - ',Outlet.Qf,'/(24*3600)*',Outlet.BSW));
        
        %Water salt mass balance
            AddEquation(obj,strcat(Inlet.ws,' - ',Outlet.ws));

		%Energy balance: isenthalpic expansion
            AddEquation(obj,strcat(...
                Inlet.Qf,'/(24*3600)*((1-',Inlet.BSW,')*',Inlet.d60,'*PP.rho_H2Ostd*',Inlet.eO, ...
						' + ',Inlet.BSW,'*',Inlet.rho_Wstd,'*',Inlet.eW,')',' - (', ...
                Outlet.Qf,'/(24*3600)*((1-',Outlet.BSW,')*',Outlet.d60,'*PP.rho_H2Ostd*',Outlet.eO, ...
						' + ',Outlet.BSW,'*',Outlet.rho_Wstd,'*',Outlet.eW,')', ...
                    ' + ',Outlet.Qg,'/(24*3600)*',Outlet.dg,'*PP.rho_AIRstd*',Outlet.eG, ...
                ')'));		

        %Pressure drop
            AddEquation(obj,strcat(DP,' - (',Inlet.P,' - ',Outlet.P,')'));

        %Pressure drop at choked flow
            AddEquation(obj,strcat(DPchok,' - ',FLpar,'^2*(',Inlet.P,' - ',FF,'*',PvL,')'));

        %Liquid critical pressure ratio factor
            AddEquation(obj,strcat(...
				FF,' - (', ...
					'0.5*(1+tanh((',PvL,'-',PpcL,')/1e2/1e-6))*(0.96 - 0.28) + ', ...
					'0.5*(1+tanh((',PpcL,'-',PvL,')/1e2/1e-6))*(0.96 - 0.28*sqrt(',PvL,'/',PpcL,'))', ...
				')'));

        %Liquid phase vapor pressure
            AddEquation(obj,strcat(PvL,' - ',Inlet.P));

        %Oil phase pseudo-critical pressure - RIAZI (2005) eq. (2.68) pp. 60
            AddEquation(obj,strcat(...
				PpcO,'/1e2',' - (', ...
					'6.9575*(',Inlet.Tbn,'^0.6791)/(',Inlet.d60,'^0.6807)', ...
					'*exp(-1.35e-2*',Inlet.Tbn,' - 0.3129*',Inlet.d60,' + 9.174e-3*',Inlet.Tbn,'*',Inlet.d60,')', ...
				')'));

        %Oil phase molar mass - API-TDB Procedure 2B2.1 (Daubert and Danner 1997)
            AddEquation(obj,strcat(...
				MwO,' - (', ...
					'20.486*(',Inlet.Tbn,'*1.8)^1.26007*(',Inlet.d60,'^4.98308)', ...
					'*exp(1.165e-4*(',Inlet.Tbn,'*1.8) - 7.78712*',Inlet.d60,' + 1.1582e-3*(',Inlet.Tbn,'*1.8)*',Inlet.d60,')', ...
				')'));

        %Molar mixing rule for liquid phase pseudo-critical pressure
            AddEquation(obj,strcat( ...
				'((1-',Inlet.BSW,')*',Inlet.d60,'*PP.rho_H2Ostd/',MwO,' + ',Inlet.BSW,'*',Inlet.rho_Wstd,'/18.02',')*',PpcL,' - (', ...
					PpcO,'*(1-',Inlet.BSW,')*',Inlet.d60,'*PP.rho_H2Ostd/',MwO,' + 220.64*1e2*',Inlet.BSW,'*',Inlet.rho_Wstd,'/18.02', ...
				')'));

        %Linear inherent characteristic function
            AddEquation(obj,strcat(fc,' - ',c));

        %Solution gas-oil ratio at equilibrium
            AddEquation(obj,strcat(...
                Rseq,' - SolGasOilRatio(PP,',Outlet.P,',',Outlet.T,',',Outlet.d60,',',Outlet.dg,')'));

        %Vaporized gas-oil ratio
            AddEquation(obj,strcat(...
                Outlet.Rv,' - VapOilGasRatio(PP,',Outlet.P,',',Outlet.T,',',Outlet.d60,',',Outlet.dg,',',Outlet.Ycg,')'));

        %Equates of outlet stream compositional variables
            AddEquation(obj,strcat(Outlet.d60,' - ',Inlet.d60));
            AddEquation(obj,strcat(Outlet.Tbn,' - ',Inlet.Tbn));
			AddEquation(obj,strcat(Outlet.dg,' - ',Inlet.dgs));
			AddEquation(obj,strcat(Outlet.Ycg,' - ',Inlet.Ycgs));

        end
    end
    
end

