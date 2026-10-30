classdef initial_head < handle
    properties
        InitVars = {};
        InitValues = [];
        NumOfInitCond = 0;
    end
    
    methods
        function [] = AddInitCondition(obj,VarName,InitValue)
            SetInitCondition(obj,VarName,InitValue);
        end
        
        function [] = SetInitCondition(obj,VarName,InitValue)
            k = GetNumOfInitCond(obj);
            SetNumOfInitCond(obj,k+1);
            obj.InitValues(k+1,1) = InitValue;
            obj.InitVars{k+1,1} = VarName;
        end
        
        function obj = SetNumOfInitCond(obj,value)
           obj.NumOfInitCond = value;
        end
        
        function value = GetNumOfInitCond(obj)
           value = obj.NumOfInitCond;
        end
    end

end

