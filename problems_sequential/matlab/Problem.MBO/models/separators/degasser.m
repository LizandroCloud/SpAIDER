classdef degasser < ModelTypes  
    properties
        L = 8;      %Vessel length          [m]
        D = 2.5;    %Vessel diameter        [m]
        V = [];     %Vessel total volume    [m^3]
    end
    
    methods
        function obj = degasser(Name,varargin)
        %MODEL BRIEF
            obj.Description = 'Horizontal vessel for gas-liquid separation';
            obj.Name = Name;

        %PARAMETERS
        %User-defined parameters
            if length(varargin) > 4
                error('Too many input parameters as arguments. Valid parameters: {L,D}');
            elseif length(varargin) > 1
                k=1;
                while k<length(varargin)
                    if strcmp(varargin{k},'L')
                        obj.L = varargin{k+1};
                    elseif strcmp(varargin{k},'D')
                        obj.D = varargin{k+1};
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
            
        %EQUIPMENT
            AddEquipment(obj,str3ph(strcat(obj.Name,'.Inlet')));
            AddEquipment(obj,liquid_stream(strcat(obj.Name,'.OutletL')));
            AddEquipment(obj,gas_stream(strcat(obj.Name,'.OutletG')));

        %VARIABLES
            AddVariable(obj,'VOLstd',volume);
            AddVariable(obj,'VOGstd',volume);
            AddVariable(obj,'VGLstd',volume);
            AddVariable(obj,'VWLstd',volume);
            AddVariable(obj,'Mws',mass);
            AddVariable(obj,'VLstd',volume);
            AddVariable(obj,'VGstd',volume);
            AddVariable(obj,'E',energy);
            AddVariable(obj,'Level',distance);

            idx = GetVarIndex(obj,'VOLstd');    obj.Variables{idx}.Name = strcat(Name,'.VOLstd');   obj.Variables{idx}.Description = 'Oil std. vol. hold-up at liquid phase';   obj.Variables{idx}.Value = 25;  obj.Variables{idx}.Unit = 'Sm^3';
            idx = GetVarIndex(obj,'VOGstd');    obj.Variables{idx}.Name = strcat(Name,'.VOGstd');   obj.Variables{idx}.Description = 'Oil std. vol. hold-up at gas phase';      obj.Variables{idx}.Value = 5;   obj.Variables{idx}.Unit = 'Sm^3';
            idx = GetVarIndex(obj,'VGLstd');    obj.Variables{idx}.Name = strcat(Name,'.VGLstd');   obj.Variables{idx}.Description = 'Gas std. vol. hold-up at liquid phase';   obj.Variables{idx}.Value = 1e3; obj.Variables{idx}.Unit = 'Sm^3';
            idx = GetVarIndex(obj,'VWLstd');    obj.Variables{idx}.Name = strcat(Name,'.VWLstd');   obj.Variables{idx}.Description = 'Water std. vol. hold-up at liquid phase'; obj.Variables{idx}.Value = 25;  obj.Variables{idx}.Unit = 'Sm^3';
            idx = GetVarIndex(obj,'Mws');       obj.Variables{idx}.Name = strcat(Name,'.Mws');      obj.Variables{idx}.Description = 'Total water salt mass hold-up';           obj.Variables{idx}.Value = 1e5; obj.Variables{idx}.Unit = 'kg';
            idx = GetVarIndex(obj,'VLstd');     obj.Variables{idx}.Name = strcat(Name,'.VLstd');    obj.Variables{idx}.Description = 'Liquid phase std. vol. hold-up';          obj.Variables{idx}.Value = 25;  obj.Variables{idx}.Unit = 'Sm^3';
            idx = GetVarIndex(obj,'VGstd');     obj.Variables{idx}.Name = strcat(Name,'.VGstd');    obj.Variables{idx}.Description = 'Gas phase std. vol. hold-up';             obj.Variables{idx}.Value = 5e3; obj.Variables{idx}.Unit = 'Sm^3';
            idx = GetVarIndex(obj,'E');         obj.Variables{idx}.Name = strcat(Name,'.E');        obj.Variables{idx}.Description = 'Energy hold-up';                          obj.Variables{idx}.Value = 5e6; obj.Variables{idx}.Unit = 'kJ';
            idx = GetVarIndex(obj,'Level');     obj.Variables{idx}.Name = strcat(Name,'.Level');    obj.Variables{idx}.Description = 'Liquid level';                            obj.Variables{idx}.Value = 1.2; obj.Variables{idx}.Unit = 'm';

        %Variable tokens
            VOLstd              = strcat(obj.Name,'.VOLstd');
            VOGstd              = strcat(obj.Name,'.VOGstd');
            VGLstd              = strcat(obj.Name,'.VGLstd');
            VWLstd              = strcat(obj.Name,'.VWLstd');
            Mws                 = strcat(obj.Name,'.Mws');
            VLstd               = strcat(obj.Name,'.VLstd');
            VGstd               = strcat(obj.Name,'.VGstd');
            Level               = strcat(obj.Name,'.Level');
            E                   = strcat(obj.Name,'.E');

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

            OutletG.Qg          = strcat(obj.Name,'.OutletG.Qg');
            OutletG.P           = strcat(obj.Name,'.OutletG.P');
            OutletG.T           = strcat(obj.Name,'.OutletG.T');
            OutletG.eG          = strcat(obj.Name,'.OutletG.eG');
            OutletG.Rv          = strcat(obj.Name,'.OutletG.Rv');
            OutletG.Bg          = strcat(obj.Name,'.OutletG.Bg');
            OutletG.d60G        = strcat(obj.Name,'.OutletG.d60G');
            OutletG.TbnG        = strcat(obj.Name,'.OutletG.TbnG');
            OutletG.dg          = strcat(obj.Name,'.OutletG.dg');
            OutletG.Ycg         = strcat(obj.Name,'.OutletG.Ycg');            
            
        %REVERSE MODEL INHERITANCE
            SystemVariables(obj);
            SystemEquations(obj);
            SystemInitConditions(obj);
        
        %EQUATIONS
        %Oil std. vol. balance
            AddEquation(obj,strcat(...
                'diff(',VOLstd,') + diff(',VOGstd,')',' - (', ...
                	Inlet.Qf,'/(24*3600)*(1-',Inlet.BSW,') + ',Inlet.Qg,'/(24*3600)*',Inlet.Rv,'*1e-3', ...
					' - ',OutletL.Qf,'/(24*3600)*(1-',OutletL.BSW,') - ',OutletG.Qg,'/(24*3600)*',OutletG.Rv,'*1e-3', ...
                ')'));

        %Oil std. vol. hold-up at liquid phase
            AddEquation(obj,strcat(VOLstd,' - ',VLstd,'*(1-',OutletL.BSW,')'));

        %Oil std. vol. hold-up at gas phase
            AddEquation(obj,strcat(VOGstd,' - ',VGstd,'*',OutletG.Rv,'*1e-3'));
            
        %Gas std. vol. balance
            AddEquation(obj,strcat(...
                'diff(',VGLstd,') + diff(',VGstd,')',' - (', ...
                	Inlet.Qf,'/(24*3600)*(1-',Inlet.BSW,')*',Inlet.Rs,' + ',Inlet.Qg,'/(24*3600)', ...
					' - ',OutletL.Qf,'/(24*3600)*(1-',OutletL.BSW,')*',OutletL.Rs,' - ',OutletG.Qg,'/(24*3600)',...
                ')'));

        %Gas std. vol. hold-up at liquid phase
            AddEquation(obj,strcat(VGLstd,' - ',VLstd,'*(1-',OutletL.BSW,')*',OutletL.Rs));
            
        %Water std. vol. balance
            AddEquation(obj,strcat(...
                'diff(',VWLstd,')',' - (',Inlet.Qf,'/(24*3600)*',Inlet.BSW,' - ',OutletL.Qf,'/(24*3600)*',OutletL.BSW,')'));
            
        %Water std. vol. hold-up
            AddEquation(obj,strcat(VWLstd,' - ',VLstd,'*',OutletL.BSW));
        
        %Water salt mass balance
            AddEquation(obj,strcat(...
                'diff(',Mws,')',' - (', ...
                Inlet.Qf,'/(24*3600)*',Inlet.BSW,'*',Inlet.rho_Wstd,'*',Inlet.ws, ...
                    ' - ',OutletL.Qf,'/(24*3600)*',OutletL.BSW,'*',OutletL.rho_Wstd,'*',OutletL.ws,...
                ')'));
        
        %Water salt mass hold-up
            AddEquation(obj,strcat(Mws,' - ',VWLstd,'*',OutletL.rho_Wstd,'*',OutletL.ws));
            
        %Volume constraint
            AddEquation(obj,strcat(...
                Vpar,' - (',VGstd,'*',OutletG.Bg,' + ',VOLstd,'*',OutletL.Bo,' + ',VWLstd,'*',OutletL.Bw,')'));

        %Liquid volume hold-up
            AddEquation(obj,strcat(...
                '(',VOLstd,'*',OutletL.Bo,' + ',VWLstd,'*',OutletL.Bw,')',' - (',...
                '(0.25*',Lpar,'*',Dpar,'^2)*(','acos(1-2*',Level,'/',Dpar,')',...
                    '- (1-2*',Level,'/',Dpar,')*sin(acos(1-2*',Level,'/',Dpar,')))',...
                ')'));

        %Energy balance
            AddEquation(obj,strcat(...
                'diff(',E,')',' - (', ...
                Inlet.Qf,'/(24*3600)*((1-',Inlet.BSW,')*',Inlet.d60,'*PP.rho_H2Ostd*',Inlet.eO, ...
                    ' + ',Inlet.BSW,'*',Inlet.rho_Wstd,'*',Inlet.eW,')', ...
                    ' + ',Inlet.Qg,'/(24*3600)*',Inlet.dg,'*PP.rho_AIRstd*',Inlet.eG,...
                ' - ',OutletL.Qf,'/(24*3600)*((1-',OutletL.BSW,')*',OutletL.d60,'*PP.rho_H2Ostd*',OutletL.eO, ...
                    ' + ',OutletL.BSW,'*',OutletL.rho_Wstd,'*',OutletL.eW,')', ...
                ' - ',OutletG.Qg,'/(24*3600)*',OutletG.dg,'*PP.rho_AIRstd*',OutletG.eG,...
                ')'));

        %Energy hold-up
            AddEquation(obj,strcat(...
                E,' - (',VOLstd,'*',OutletL.d60,'*PP.rho_H2Ostd*',OutletL.eO, ...
                    ' + ',VWLstd,'*',OutletL.rho_Wstd,'*',OutletL.eW,' + ',VGstd,'*',OutletG.dg,'*PP.rho_AIRstd*',OutletG.eG,' - ',OutletG.P,'*',Vpar,...
                ')'));
            
        %Property mixing rule for oil specific gravity (std. vol. balance)
            AddEquation(obj,strcat(...
                OutletL.d60,'*diff(',VOLstd,') + ',VOLstd,'*diff(',OutletL.d60,') + ',OutletG.d60G,'*diff(',VOGstd,') + ',VOGstd,'*diff(',OutletG.d60G,')',' - (', ...
                Inlet.Qf,'/(24*3600)*(1-',Inlet.BSW,')*',Inlet.d60,' + ',Inlet.Qg,'/(24*3600)*',Inlet.Rv,'*1e-3','*',Inlet.d60G,...
                    ' - ',OutletL.Qf,'/(24*3600)*(1-',OutletL.BSW,')*',OutletL.d60,' - ',OutletG.Qg,'/(24*3600)*',OutletG.Rv,'*1e-3','*',OutletG.d60G, ...
                ')'));

        %Property mixing rule for Tbn by means of Watson K-factor (std. vol. balance)
            AddEquation(obj,strcat(...
                '(',OutletL.Tbn,')^(1/3)','*diff(',VOLstd,') + ',VOLstd,'*1/3*',OutletL.Tbn,'^(-2/3)*diff(',OutletL.Tbn,') + ',...
                    '(',OutletG.TbnG,')^(1/3)','*diff(',VOGstd,') + ',VOGstd,'*1/3*',OutletG.TbnG,'^(-2/3)*diff(',OutletG.TbnG,')',' - (', ...
                	Inlet.Qf,'/(24*3600)*(1-',Inlet.BSW,')*(',Inlet.Tbn,')^(1/3)', ...
					' + ',Inlet.Qg,'/(24*3600)*',Inlet.Rv,'*1e-3','*(',Inlet.TbnG,')^(1/3)',...
                    ' - ',OutletL.Qf,'/(24*3600)*(1-',OutletL.BSW,')*(',OutletL.Tbn,')^(1/3)',...
                    ' - ',OutletG.Qg,'/(24*3600)*',OutletG.Rv,'*1e-3','*(',OutletG.TbnG,')^(1/3)', ...
                ')'));

        %Property mixing rule for gas specific gravity (std vol. balance)
            AddEquation(obj,strcat(...
                OutletL.dg,'*diff(',VGLstd,') + ',VGLstd,'*diff(',OutletL.dg,') + ', ...
					OutletG.dg,'*diff(',VGstd,') + ',VGstd,'*diff(',OutletG.dg,')',' - (', ...
                	Inlet.Qf,'/(24*3600)*(1-',Inlet.BSW,')*',Inlet.Rs,'*',Inlet.dgs,' + ',Inlet.Qg,'/(24*3600)*',Inlet.dg,...
                    ' - ',OutletL.Qf,'/(24*3600)*(1-',OutletL.BSW,')*',OutletL.Rs,'*',OutletL.dg,' - ',OutletG.Qg,'/(24*3600)*',OutletG.dg,...
                ')'));

        %Property mixing rule for condensate/gas yield (std. vol. balance)
            AddEquation(obj,strcat(...
                OutletL.Ycg,'*diff(',VGLstd,') + ',VGLstd,'*diff(',OutletL.Ycg,') + ', ...
					OutletG.Ycg,'*diff(',VGstd,') + ',VGstd,'*diff(',OutletG.Ycg,')',' - (', ...
                	Inlet.Qf,'/(24*3600)*(1-',Inlet.BSW,')*',Inlet.Rs,'*',Inlet.Ycgs,' + ',Inlet.Qg,'/(24*3600)*',Inlet.Ycg,...
                    ' - ',OutletL.Qf,'/(24*3600)*(1-',OutletL.BSW,')*',OutletL.Rs,'*',OutletL.Ycg,' - ',OutletG.Qg,'/(24*3600)*',OutletG.Ycg,...
                ')'));
            
        %Solution gas-oil ratio
            AddEquation(obj,strcat(...
                OutletL.Rs,' - SolGasOilRatio(PP,',OutletG.P,',',OutletL.T,',',OutletL.d60,',',OutletG.dg,')'));

        %Vaporized gas-oil ratio
            AddEquation(obj,strcat(...
                OutletG.Rv,' - VapOilGasRatio(PP,',OutletG.P,',',OutletG.T,',',OutletL.d60,',',OutletG.dg,',',OutletG.Ycg,')'));

        %Hypothetical oil specific gravity at gas phase
            AddEquation(obj,strcat(...
                OutletG.d60G,' - HypOilSpecGrav(PP,',OutletG.P,',',OutletG.T,',',OutletL.d60,')'));

        %Hypothetical oil normal boiling point at gas phase
            AddEquation(obj,strcat(...
                OutletG.TbnG,' - HypOilNormBolPnt(PP,',OutletG.P,',',OutletG.T,',',OutletL.d60,',',OutletL.Tbn,')'));

        %Solution gas specific gravity
            AddEquation(obj,strcat(...
                OutletL.dg,' - SolGasSpecGrav(PP,',OutletG.P,',',OutletG.T,',',OutletG.dg,')'));
            
        %Solution condensate-gas yield
            AddEquation(obj,strcat(...
                OutletL.Ycg,' - SolCondGasRatio(PP,',OutletG.P,',',OutletG.T,',',OutletG.Ycg,')'));

        %Thermal equlibrium
            AddEquation(obj,strcat(OutletL.T,' - ',OutletG.T));

        %Mechanical equlibrium
            AddEquation(obj,strcat(Inlet.P,' - ',OutletG.P));
            AddEquation(obj,strcat(OutletL.P,' - ',OutletG.P));

        end
    end
    
end

