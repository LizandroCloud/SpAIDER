classdef triphase_stream < str3ph
    methods
        function obj = triphase_stream(Name)
       %MODEL BRIEF
            obj = obj@str3ph(Name);
            obj.Description = '3-phase stream model for physical property calculations';
            obj.Name = Name;

       %VARIABLES
            AddVariable(obj,'Zg',admensional);
            
            idx = GetVarIndex(obj,'Zg');    obj.Variables{idx}.Name = strcat(Name,'.Zg');   obj.Variables{idx}.Description = 'Gas compressibility factor'; obj.Variables{idx}.Value = 1;    obj.Variables{idx}.Unit = '';
            
       %EQUATIONS
       %Variable tokens
            BSW = strcat(obj.Name,'.BSW');
            P = strcat(obj.Name,'.P');
            T = strcat(obj.Name,'.T');
            eO = strcat(obj.Name,'.eO');
            eG = strcat(obj.Name,'.eG');
            eW = strcat(obj.Name,'.eW');
            rhoL = strcat(obj.Name,'.rhoL');
            rhoG = strcat(obj.Name,'.rhoG');
            rho_Wstd = strcat(obj.Name,'.rho_Wstd');
            Bo = strcat(obj.Name,'.Bo');
            Bg = strcat(obj.Name,'.Bg');
            Bw = strcat(obj.Name,'.Bw');
            d60 = strcat(obj.Name,'.d60');
            d60G = strcat(obj.Name,'.d60G');
            Tbn = strcat(obj.Name,'.Tbn');
            TbnG = strcat(obj.Name,'.TbnG');
            dg = strcat(obj.Name,'.dg');
            dgs = strcat(obj.Name,'.dgs');
            Ycg = strcat(obj.Name,'.Ycg');
            Ycgs = strcat(obj.Name,'.Ycgs');
            ws = strcat(obj.Name,'.ws');

            Zg = strcat(obj.Name,'.Zg');
            
        %Oil specific enthalpy
            AddEquation(obj,strcat(...
                eO,' - OilSpecEnth(PP,',T,',',d60,',',Tbn,')'));

        %Gas specific enthalpy
            AddEquation(obj,strcat(...
                eG,' - (IdealGasSpecEnth(PP,',T,',',dg,') + GasLatHeat(PP,',P,',',T,',',d60,',',dg,') + GasSpecResEnth(PP,',P,',',T,',',dg,'))'));
        
        %Water specific enthalpy    
            AddEquation(obj,strcat(...
                eW,' - WaterSpecEnth(PP,',P,',',T,',',ws,')'));

        %Liquid phase mass density
            AddEquation(obj,strcat(...
                rhoL,'*((1-',BSW,')*',Bo,'+ ',BSW,'*',Bw,') - ((1-',BSW,')*',d60,'*PP.rho_H2Ostd + ',BSW,'*',rho_Wstd,')'));
        
        %Gas phase mass density
            AddEquation(obj,strcat(...
                rhoG,'*',Zg,'*PP.R*',T,'-(',P,'*PP.MM_Air*',dg,')'));

        %Water mass density at std. cond.
            AddEquation(obj,strcat(...
                rho_Wstd,'- WaterStdMassDens(PP,',ws,')'));
            
        %Oil phase formation volume factor
            AddEquation(obj,strcat(...
                Bo,'- OilFVF(PP,',P,',',T,',',d60,',',dg,')'));
            
        %Gas phase formation volume factor
            AddEquation(obj,strcat(...
                Bg,'-(PP.Pref/',P,')*(',T,'/PP.Tref)*',Zg));
            
        %Gas compressibility factor
            AddEquation(obj,strcat(...
                Zg,'- EvalGasCompFactor(PP,',P,',',T,',',dg,',',Zg,')'));

        %Water phase formation volume factor
            AddEquation(obj,strcat(...
                Bw,'- WaterFVF(PP,',P,',',T,')'));
            
        %Hypothetical oil specific gravity at gas phase
            AddEquation(obj,strcat(...
                d60G,'- HypOilSpecGrav(PP,',P,',',T,',',d60,')'));

        %Hypothetical oil normal boiling point at gas phase
            AddEquation(obj,strcat(...
                TbnG,'- HypOilNormBolPnt(PP,',P,',',T,',',d60,',',Tbn,')'));

        %Solution gas specific gravity
            AddEquation(obj,strcat(...
                dgs,'- SolGasSpecGrav(PP,',P,',',T,',',dg,')'));
            
        %Solution condensate-gas yield
            AddEquation(obj,strcat(...
                Ycgs,'- SolCondGasRatio(PP,',P,',',T,',',Ycg,')'));
            
        end
    end
end