classdef liquid_stream < strL
    methods
        function obj = liquid_stream(Name)
        %MODEL BRIEF
            obj = obj@strL(Name);
            obj.Description = 'Liquid (oil + water) stream model for physical property calculations';
            obj.Name = Name;

       %EQUATIONS
       %Variable tokens
            BSW = strcat(obj.Name,'.BSW');
            P = strcat(obj.Name,'.P');
            T = strcat(obj.Name,'.T');
            eO = strcat(obj.Name,'.eO');
            eW = strcat(obj.Name,'.eW');
            rhoL = strcat(obj.Name,'.rhoL');
            rho_Wstd = strcat(obj.Name,'.rho_Wstd');
            Bo = strcat(obj.Name,'.Bo');
            Bw = strcat(obj.Name,'.Bw');
            d60 = strcat(obj.Name,'.d60');
            Tbn = strcat(obj.Name,'.Tbn');
            dg = strcat(obj.Name,'.dg');
            ws = strcat(obj.Name,'.ws');
            
        %Oil specific enthalpy
            AddEquation(obj,strcat(...
                eO,' - OilSpecEnth(PP,',T,',',d60,',',Tbn,')'));
        
        %Water specific enthalpy    
            AddEquation(obj,strcat(...
                eW,' - WaterSpecEnth(PP,',P,',',T,',',ws,')'));

        %Liquid phase mass density
            AddEquation(obj,strcat(...
                rhoL,'*((1-',BSW,')*',Bo,'+ ',BSW,'*',Bw,') - ((1-',BSW,')*',d60,'*PP.rho_H2Ostd + ',BSW,'*',rho_Wstd,')'));

        %Water mass density at std. cond.
            AddEquation(obj,strcat(...
                rho_Wstd,'- WaterStdMassDens(PP,',ws,')'));
            
        %Oil phase formation volume factor
            AddEquation(obj,strcat(...
                Bo,'- OilFVF(PP,',P,',',T,',',d60,',',dg,')'));

        %Water phase formation volume factor
            AddEquation(obj,strcat(...
                Bw,'- WaterFVF(PP,',P,',',T,')'));
            
        end
    end
end

