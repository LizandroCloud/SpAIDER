classdef strG < strGL
    methods
        function obj = strG(Name)
        %MODEL BRIEF
            obj = obj@strGL(Name);
            obj.Description = 'Gas only stream basic model';
            obj.Name = Name;

        %EQUATIONS
        %Variables tokens
            Qo = strcat(obj.Name,'.Qo');
            eO = strcat(obj.Name,'.eO');
            rhoL = strcat(obj.Name,'.rhoL');
            Rs = strcat(obj.Name,'.Rs');
            Bo = strcat(obj.Name,'.Bo');
            d60 = strcat(obj.Name,'.d60');
            d60G = strcat(obj.Name,'.d60G');
            Tbn = strcat(obj.Name,'.Tbn');
            TbnG = strcat(obj.Name,'.TbnG');
            dg = strcat(obj.Name,'.dg');
            dgs = strcat(obj.Name,'.dgs');
            Ycg = strcat(obj.Name,'.Ycg');
            Ycgs = strcat(obj.Name,'.Ycgs');
        
        %No liquid: oil std. vol. flow rate
            AddEquation(obj,strcat(Qo,'-0.00'));
        
        %No liquid: oil specific enthalpy
            AddEquation(obj,strcat(eO,'-0.00'));

        %No liquid: liquid phase mass density
            AddEquation(obj,strcat(rhoL,'-1000'));

        %No liquid: solution gas-oil ratio
            AddEquation(obj,strcat(Rs,'-0.00'));

        %No liquid: oil phase formation volume factor
            AddEquation(obj,strcat(Bo,'-1.00'));
            
        %No liquid: hypothetical oil specific gravity 60 F
            AddEquation(obj,strcat(d60G,'-',d60));

        %No liquid: hypothetical oil normal boiling point
            AddEquation(obj,strcat(TbnG,'-',Tbn));

        %No liquid: solution gas specific gravity
            AddEquation(obj,strcat(dgs,'-',dg));

        %No liquid: condensate-gas yield
            AddEquation(obj,strcat(Ycgs,'-',Ycg));
        end
    end
end

