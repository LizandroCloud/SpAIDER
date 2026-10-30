classdef dens_mass < VariableTypes    
    methods
        function obj = dens_mass(~)
            obj.Value = 1000;
            obj.Description = 'Mass density';
            obj.Unit = 'kg/m^3';
        end
    end    
end

