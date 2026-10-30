classdef variables_head < handle
    properties
        Variables = {};
        NumOfVars = 0;
    end
    
    methods
        function [] = AddVariable(obj,VarName,VarType)
            SetVariable(obj,VarName,VarType);
        end
        
        function [] = SetVariable(obj,VarName,VarType)
            k = GetNumOfVars(obj);
            SetNumOfVars(obj,k+1);
            obj.Variables{k+1,1} = VarType;
            obj.Variables{k+1}.Name = VarName;
%            obj.Variables{k+1,1}.ID = k+1;
        end
        
        function obj = SetNumOfVars(obj,value)
           obj.NumOfVars = value;
        end
        
        function value = GetNumOfVars(obj)
           value = obj.NumOfVars;
        end
        
        function idx = GetVarIndex(obj,VarName)
            torf = true;
            for k=1:length(obj.Variables)
                if strcmp(obj.Variables{k}.Name,VarName) == true && torf
                    idx = k;
                    torf = false;
                end
            end
        end
        
    end

end

