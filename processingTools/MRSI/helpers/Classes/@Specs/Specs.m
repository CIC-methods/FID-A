classdef Specs < Voxel
    % Child class of Voxel.Contains spectrum information as well
    
    properties
        specs
    end
    
    methods
        function obj = Voxel(x, y, z, index)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj = Voxel(x, y, z, index);
        end
        
        function obj = set_specs(obj,specs)
            %set_specs. Sets the specs property of the object
            obj.specs = specs;
        end
    end
end

