classdef volume < VariableTypes    
    methods
        function obj = volume(~)
            obj.Value = 1e3;
            obj.Description = 'Volume';
            obj.Unit = 'm^3';
        end
    end    
end

