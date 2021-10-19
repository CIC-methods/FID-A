%op_CSIDensityCompensation.m
%Brenden Kadota, Jamie Near, Sunnybrook 2021.
%
% USAGE:
% [out, weight_map]=op_CSIDensityCompensation(in,k_space_file, NamveValueArgs);
% 
% DESCRIPTION:
% Takes a CSI FID-A structure and a k space file and applies density
% compensation from voronoi diagrams. Name value arguments can be passed to
% the function to plot the voronoi diagrams or chand the area of support.
% 
% INPUTS:
% in        = input data in FID-A CSI structure format.
% k_space_file     = k space file with the first two columns being kx and ky. 
% NAME VALUE ARGUMENTS:
%      Name: 'isPlot' Values: true of false. Plots vornoi diagram
%      Name: 'areaOfSupport' Values: 'circular', 'rectangular'. Adds points
%      around k trajectory in the shape selected.


function [out, weight_matrix] = op_CSIDensityCompensation(in, k_space_file, NameValueArgs)
    arguments
        in (1,1) struct
        k_space_file (1,:) char {mustBeFile}
        NameValueArgs.isPlot (1,1) logical = false
        NameValueArgs.areaOfSupport (1, :) {mustBeMember(NameValueArgs.areaOfSupport,{'rectangular','circular'})} = 'circular'
    end
    
    k_table = readtable(k_space_file);
    k_space = [k_table.K_x, k_table.K_y];
    
    [index, v_outer] = convhull(k_space(:,1), k_space(:,2));
    inner_traj = k_space;
    inner_traj(index, :) = [];
    [~, v_inner] = convhull(inner_traj(:,1), inner_traj(:,2));
    
    alpha = sqrt(v_outer/v_inner);
    %get max diameter of k space
    diameter = max(k_space,[] , 'all') - min(k_space,[] , 'all');

    if(NameValueArgs.areaOfSupport == "circular")
        %equal spaced points around 2pi
        theta = 0:0.01:2*pi-0.01;
        %create  a border around k_space in complex plane
        border = (diameter*alpha/2)*exp(1i*theta);
        border = [real(border)' imag(border)'];
    elseif(NameValueArgs.areaOfSupport == "rectangular")
        width = diameter*alpha;
        edge = linspace(-width/2,width/2,250);
        top = [edge' repmat(width/2, size(edge,2),1)];
        right = [repmat(width/2, size(edge,2),1), edge'];
        left = [repmat(-width/2, size(edge,2),1), edge'];
        bottom = [edge' repmat(-width/2, size(edge,2),1)];
        border = cat(1, top, right,left, bottom);
    end
    %concat border with k_space
    bordered_k_space = cat(1, k_space, border);
    
    [unique_k,~,IC] = unique(bordered_k_space, 'rows', 'stable');
    [occurences, indecies] = groupcounts(IC);
    %get delaunay triangulation from points
    tr = delaunayTriangulation(unique_k(:,1), unique_k(:,2));
    
    %create voronoi diagram from triangulation
    [vertecies,edge_indecies] = voronoiDiagram(tr);
    %preallocation
    vol = zeros(size(edge_indecies,1), 1);
    
    %loop through 
    for i = 1:length(edge_indecies)
        %This happends when voronoi diagram goes to inf
        if(all(edge_indecies{i}~=1))
            %Calculate the volume from Voronoi diagram 
            [~, vol(i)] = convhulln(vertecies(edge_indecies{i}, :));
        end
    end
    %apply weights for duplicate points
    final_vol = arrayfun(@(occ, ind) vol(ind)/occ , occurences, indecies);
    %remove zeros (convex hull)
    weight_points_index = find(final_vol);
    final_vol = final_vol(weight_points_index);
    final_vol = final_vol./max(final_vol, [], 'all');

    %calculate density (rho_i = 1/w_i);
    density = 1./final_vol;


    point_indecies = IC(1:size(k_space,1));
    weights = arrayfun(@(x) final_vol(IC(x)), point_indecies);
    num_TR = max(k_table.TR, [], 'all');
    spatial_points = height(k_table)/num_TR;
    weight_matrix = reshape(weights, [spatial_points, num_TR]);


    for i = 1:spatial_points:in.sz(1)
        if(i+spatial_points-1 > in.sz(1))
            in.fids(i:end,:) = in.fids(i:end,:).*weight_matrix(in.sz(1)-i+1,:);
        else 
            in.fids(i:i+spatial_points-1,:) = in.fids(i:i+spatial_points-1,:).*weight_matrix;
        end
    end

    out = in;


    if(NameValueArgs.isPlot)
        figure
        %plotting the sampling points
        subplot(2,2,1)
        scatter(tr.Points(:,1), tr.Points(:,2)), title('Sampling Points with Border');

        %plotting voronoi diagram
        subplot(2,2,2)
        hold on
        scatter(tr.Points(:,1), tr.Points(:,2), 30, '.'), title('Voronoi diagram');
        edges = cellfun(@(point) (vertecies(point,:)), edge_indecies, 'UniformOutput', false);
        cellfun(@(convex_hull) plot(convex_hull(:,1), convex_hull(:,2), '-r'), edges);
        hold off

        subplot(2,2,3)
        x = tr.Points(weight_points_index,1);
        y = tr.Points(weight_points_index,2);

        scatter(x, y, [], final_vol);
        colorbar

        subplot(2,2,4)
        x = tr.Points(weight_points_index,1);
        y = tr.Points(weight_points_index,2);
        [x_interp,y_interp] = meshgrid(-diameter/2:1/(4*in.fovX):diameter/2,-diameter/2:1/(4*in.fovY):diameter/2);
        zi = griddata(x,y,density,x_interp,y_interp);

        mesh(x_interp, y_interp, zi);
        colorbar
    end

end
