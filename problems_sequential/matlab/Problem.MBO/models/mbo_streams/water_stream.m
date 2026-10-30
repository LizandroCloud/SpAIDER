classdef water_stream < strW
    methods
        function obj = water_stream(Name)
        %MODEL BRIEF
            obj = obj@strW(Name);
            obj.Description = 'Water stream model for physical property calculations';
            obj.Name = Name;

       %EQUATIONS
       %Variable tokens
            Qw = strcat(obj.Name,'.Qw');
			TOG = strcat(obj.Name,'.TOG');
            P = strcat(obj.Name,'.P');
            T = strcat(obj.Name,'.T');
            eW = strcat(obj.Name,'.eW');
            rhoW = strcat(obj.Name,'.rhoW');
            rho_Wstd = strcat(obj.Name,'.rho_Wstd');
            d60 = strcat(obj.Name,'.d60');
            Tbn = strcat(obj.Name,'.Tbn');
            ws = strcat(obj.Name,'.ws');
        
        %Water specific enthalpy    
            AddEquation(obj,strcat(...
                eW,' - WaterSpecEnth(PP,',P,',',T,',',ws,')'));

        %Liquid phase mass density
            AddEquation(obj,strcat(...
                rhoW,' - WaterFVF(PP,',P,',',T,')*',rho_Wstd));

        %Water mass density at std. cond.
            AddEquation(obj,strcat(...
                rho_Wstd,'- WaterStdMassDens(PP,',ws,')'));
        end
    end
end

