classdef strL < str3ph
    methods
        function obj = strL(Name)
        %MODEL BRIEF
            obj = obj@str3ph(Name);
            obj.Description = 'Liquid (oil + water) stream basic model';
            obj.Name = Name;

        %EQUATIONS
        %Variables tokens
            Qg = strcat(obj.Name,'.Qg');
            eG = strcat(obj.Name,'.eG');
            rhoG = strcat(obj.Name,'.rhoG');
            Rv = strcat(obj.Name,'.Rv');
            Bg = strcat(obj.Name,'.Bg');
            d60 = strcat(obj.Name,'.d60');
            d60G = strcat(obj.Name,'.d60G');
            Tbn = strcat(obj.Name,'.Tbn');
            TbnG = strcat(obj.Name,'.TbnG');
            dg = strcat(obj.Name,'.dg');
            dgs = strcat(obj.Name,'.dgs');
            Ycg = strcat(obj.Name,'.Ycg');
            Ycgs = strcat(obj.Name,'.Ycgs');
        
        %No gas: gas std. vol. flow rate
            AddEquation(obj,strcat(Qg,'-0.00'));
        
        %No gas: gas specific enthalpy
            AddEquation(obj,strcat(eG,'-0.00'));

        %No gas: gas phase mass density
            AddEquation(obj,strcat(rhoG,'-1.00'));

        %No gas: vaporized oil-gas ratio
            AddEquation(obj,strcat(Rv,'-0.00'));
            
        %No gas: gas phase formation volume factor
            AddEquation(obj,strcat(Bg,'-1.00'));

        %No gas: hypothetical oil specific gravity 60 F
            AddEquation(obj,strcat(d60G,'-',d60));

        %No gas: hypothetical oil normal boiling point
            AddEquation(obj,strcat(TbnG,'-',Tbn));

        %No gas: solution gas specific gravity
            AddEquation(obj,strcat(dgs,'-',dg));

        %No gas: condensate-gas yield
            AddEquation(obj,strcat(Ycgs,'-',Ycg));
        end
    end
end