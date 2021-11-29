classdef Voxel 
    properties
        %points is a [3, 8] matrix with the first dimensions representign
        %x,y,z coordinates and the second dimension representing
        %individual points
        points
        %Edges connects points by indexes ie (1,2) connects point 1 to
        %point 2
        edges
        %center: Center point of the voxel
        center
        %Index: [3,1] vector mapping back to position in MRSI object's FIDS
        %or SPECS properties
        index
        %minimum values along each x,y,z coordinate.(for speedup)
        min_coords
        %maximum values along each x,y,z coordinate. (for speedup)
        max_coords
    end
    methods
        function obj = Voxel(x, y, z, index)
            if(exist('x','var'))
                [mesh_x, mesh_y, mesh_z] = meshgrid([x(1), x(2)], [y(1), y(2)], [z(1),z(2)]);
                points = [mesh_x(:) mesh_y(:), mesh_z(:)];
                edges = [1,2; 1,3; 3,4; 2,4; 1,5; 2,6; 3,7; 4,8; 5,6; 5,7; 6,8; 7,8];
                obj.edges = edges;
                obj.points = points';
                obj.center = mean(obj.points, 2);
                obj.index = index;
                obj.min_coords = min(obj.points, [], 2);
                obj.max_coords = max(obj.points, [], 2);
            end
        end

        
        % is_intersect
        % returns true if the plane intersects with the voxel, returns
        % false otherwise. Plane is along cartesian axes and is defined by
        % dimension number and coordinate of that dimension. The plane is then defined along the other two axes. 
        % ie.) plane_dims = 3 (z dimension), coord = 200 (value of
        % coordinate), plane is at z = 200 along x and y axes (transverse
        % plane)
        %%%%%%%%%%%%%%%%%%
        % NOT USED ANYMORE
        %%%%%%%%%%%%%%%%%%
        %function bool = is_intersect(obj, plane_dims, coord)
        %    min_coord = obj.min_coords(plane_dims);
        %   max_coord = obj.max_coords(plane_dims);
        %    bool = coord > min_coord && coord < max_coord;
        %end
        
        %find_intersection
        %Finds the points that are intersecting with the plane. The points
        %are then given in sequential order to be plotted, with the first
        %point beign given a second time.
        function box = find_intersection(obj, plane_dims, coord)
            box = zeros(2,4);
            counter = 1;
            
            set_diff = [1,2,3];
            set_diff = set_diff(~(set_diff == plane_dims));
            for i = 1:size(obj.edges,1)
                vox_points = obj.points(:,obj.edges(i,:));
                first_point = vox_points(:,1);
                second_point = vox_points(:,2);
                smaller = min(first_point(plane_dims), second_point(plane_dims));
                bigger = max(first_point(plane_dims), second_point(plane_dims));
                if(coord > smaller && coord < bigger)
                    diff = first_point - second_point;
                    length = (coord - second_point(plane_dims))/ diff(plane_dims);
                    intersect = diff*length + second_point;
                    
                    box(:, counter) = intersect(set_diff);
                    counter = counter + 1;
                end
            end
            if(counter == 1)
                box = [];
            else
                c = mean(box, 2);
                d = box - c;
                th = atan2(d(2,:),d(1,:)); % angle above x axis
                [~, idx] = sort(th);
                box = box(:, idx);
                box = [box box(:,1)];
            end
            
        end
    end

end
