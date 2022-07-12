classdef Coordinate
    %COORDINATE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Access = private)
        x
        y
        z
    end
    
    methods
        function obj = Coordinate(x, y, z)
            %COORDINATE Construct an instance of this class
            %   Detailed explanation goes here
            obj.x = x;
            obj.y = y;
            obj.z = z;
        end
        
        function xCoordinate = getX(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            xCoordinate = obj.x;
        end

        function xCoordinate = getY(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            xCoordinate = obj.y;
        end

        function xCoordinate = getZ(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            xCoordinate = obj.z;
        end
    end
end

