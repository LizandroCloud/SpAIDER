classdef gas_stream < strG
    methods
        function obj = gas_stream(Name)
        %MODEL BRIEF
            obj = obj@strG(Name);
            obj.Description = 'Gas only stream model for physical property calculations';
            obj.Name = Name;

       %VARIABLES
            AddVariable(obj,'Zg',admensional);
            
            idx = GetVarIndex(obj,'Zg');    obj.Variables{idx}.Name = strcat(Name,'.Zg');   obj.Variables{idx}.Description = 'Gas compressibility factor'; obj.Variables{idx}.Value = 1;    obj.Variables{idx}.Unit = '';
            
       %EQUATIONS
       %Variable tokens
            P = strcat(obj.Name,'.P');
            T = strcat(obj.Name,'.T');
            eG = strcat(obj.Name,'.eG');
            rhoG = strcat(obj.Name,'.rhoG');
            Bg = strcat(obj.Name,'.Bg');
            d60 = strcat(obj.Name,'.d60');
            dg = strcat(obj.Name,'.dg');
            Zg = strcat(obj.Name,'.Zg');

        %Gas specific enthalpy
            AddEquation(obj,strcat(...
                eG,' - (IdealGasSpecEnth(PP,',T,',',dg,') + GasLatHeat(PP,',P,',',T,',',d60,',',dg,') + GasSpecResEnth(PP,',P,',',T,',',dg,'))'));
        
        %Gas phase mass density
            AddEquation(obj,strcat(...
                rhoG,'*',Zg,'*PP.R*',T,'-(',P,'*PP.MM_Air*',dg,')'));

        %Gas phase formation volume factor
            AddEquation(obj,strcat(...
                Bg,'-(PP.Pref/',P,')*(',T,'/PP.Tref)*',Zg));
            
        %Gas compressibility factor
            AddEquation(obj,strcat(...
                Zg,'- EvalGasCompFactor(PP,',P,',',T,',',dg,',',Zg,')'));
        end
    end
end

