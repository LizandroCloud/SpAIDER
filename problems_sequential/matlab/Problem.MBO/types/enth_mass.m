classdef enth_mass < VariableTypes    
    methods
        function obj = enth_mass(~)
            obj.Value = 250;
            obj.Description = 'Specific enthalpy';
            obj.Unit = 'kJ/kg';
        end
    end    
end

