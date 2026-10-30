classdef liquid_source < ModelTypes
    methods
        function obj = liquid_source(Name)
        %MODEL BRIEF
            obj.Description = 'Liquid source stream model';
            obj.Name = Name;

        %EQUIPMENT
            AddEquipment(obj,liquid_stream(strcat(obj.Name,'.Outlet')));

        %REVERSE MODEL INHERITANCE
            SystemVariables(obj);
            SystemEquations(obj);
            SystemInitConditions(obj);
            
        %EQUATIONS
       %Variable tokens
            Outlet.P = strcat(obj.Name,'.Outlet.P');
            Outlet.T = strcat(obj.Name,'.Outlet.T');
            Outlet.Rs = strcat(obj.Name,'.Outlet.Rs');
            Outlet.d60 = strcat(obj.Name,'.Outlet.d60');
            Outlet.dg = strcat(obj.Name,'.Outlet.dg');
        
        %Solution gas-oil ratio
            AddEquation(obj,strcat(...
                Outlet.Rs,'- SolGasOilRatio(PP,',Outlet.P,',',Outlet.T,',',Outlet.d60,',',Outlet.dg,')'));
        end
    end
end

