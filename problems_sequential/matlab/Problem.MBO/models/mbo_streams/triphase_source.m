classdef triphase_source < ModelTypes
    methods
        function obj = triphase_source(Name)
        %MODEL BRIEF
            obj.Description = '3-phase source stream model';
            obj.Name = Name;

        %EQUIPMENT
            AddEquipment(obj,triphase_stream(strcat(obj.Name,'.Outlet')));

        %REVERSE MODEL INHERITANCE
            SystemVariables(obj);
            SystemEquations(obj);
            SystemInitConditions(obj);
            
        %EQUATIONS
       %Variable tokens
            Outlet.P = strcat(obj.Name,'.Outlet.P');
            Outlet.T = strcat(obj.Name,'.Outlet.T');
            Outlet.Rs = strcat(obj.Name,'.Outlet.Rs');
            Outlet.Rv = strcat(obj.Name,'.Outlet.Rv');
            Outlet.d60 = strcat(obj.Name,'.Outlet.d60');
            Outlet.dg = strcat(obj.Name,'.Outlet.dg');
            Outlet.Ycg = strcat(obj.Name,'.Outlet.Ycg');
        
        %Solution gas-oil ratio
            AddEquation(obj,strcat(...
                Outlet.Rs,'- SolGasOilRatio(PP,',Outlet.P,',',Outlet.T,',',Outlet.d60,',',Outlet.dg,')'));
        
        %Vaporized gas-oil ratio
            AddEquation(obj,strcat(...
                Outlet.Rv,'- VapOilGasRatio(PP,',Outlet.P,',',Outlet.T,',',Outlet.d60,',',Outlet.dg,',',Outlet.Ycg,')'));
        end
    end
end

