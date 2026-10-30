classdef flow_vol < VariableTypes    
    methods
        function obj = flow_vol(~)
            obj.Value = 1e4;
            obj.Description = 'Volumetric flow rate';
            obj.Unit = 'm^3/d';
        end
    end    
end

