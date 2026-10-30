classdef molweight < VariableTypes    
    methods
        function obj = molweight(~)
            obj.Value = 20;
            obj.Description = 'Molar mass';
            obj.Unit = 'kg/kmol';
        end
    end    
end

