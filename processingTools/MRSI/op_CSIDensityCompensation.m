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


function [MRSIStruct, weightMatrix] = op_CSIDensityCompensation(MRSIStruct, k_space_file, NameValueArgs)
    arguments
        MRSIStruct (1,1) struct
        k_space_file (1,:) char {mustBeFile}
        NameValueArgs.isPlot (1,1) logical = false
        NameValueArgs.areaOfSupport (1, :) {mustBeMember(NameValueArgs.areaOfSupport,{'rectangular','circular'})} = 'circular'
    end
    if(getFlags(MRSIStruct, 'spatialft') == true)
        error('Please keep the image data in the k-space when using density compensation')
    end
    k_space = getKSpace(k_space_file);

    weights = voronoi(k_space_file, NameValueArgs, k_space);
    %calculate density (desnity = 1/w_i);


    spatialPoints = calculateNumSpatialPoints(k_space_file);
    timePoints = floor(getSizeFromDimensions(MRSIStruct, {'t'}) / spatialPoints);
    TR = length(k_space)/spatialPoints;
    
    

    %permute so only t and y dimensions are first
    [MRSIStruct, prev_permute, prev_size] = reshapeDimensions(MRSIStruct, {'t', 'ky'});
    data = getData(MRSIStruct);
    %get normal size to reshape back to
    for iTime = 1:timePoints
        startTime = spatialPoints * (iTime - 1) + 1;
        endTime = spatialPoints * iTime;
        data(startTime:endTime, :, :) = data(startTime:endTime, :, :) .* weights;
    end
    MRSIStruct = setData(MRSIStruct, data);
    MRSIStruct = reshapeBack(MRSIStruct, prev_permute, prev_size);


end

function kSpace = getKSpace(filename)
    k_table = readtable(filename);
    kSpace = [k_table.Kx, k_table.Ky];
end

function diameter = getMaxDiameterofKSpace(fileName)
    kSpace = getKSpace(fileName);
    diameter = max(kSpace,[] , 'all') - min(kSpace,[] , 'all');
end

function border = createBorder(fileName, areaOfSupport)
    kSpace = getKSpace(fileName);
    
    %get convex hull
    [convexHullIndex, volumeOuter] = convhull(kSpace(:,1), kSpace(:,2));
    inner_traj = kSpace;
    %remove convex hull points
    inner_traj(convexHullIndex, :) = [];
    %get new convex hull
    [~, volumeInner] = convhull(inner_traj(:,1), inner_traj(:,2));
    
    %difference of volume of both convex hulls becomes alpha
    alpha = volumeOuter/volumeInner;
    %get max diameter of k space
    
    diameter = getMaxDiameterofKSpace(fileName);
    if(strcmp(areaOfSupport, "circular"))
        %equal spaced points around 2pi
        theta = 0:0.005:2*pi-0.005;
        %create  a border around k_space in complex plane
        border = (diameter*alpha/2)*exp(1i*theta);
        border = [real(border)' imag(border)'];
    elseif(strcmp(areaOfSupport, "rectangular"))
        width = diameter*alpha;
        edge = linspace(-width/2,width/2,250);
        
        top = [edge' repmat(width/2, size(edge,2),1)];
        right = [repmat(width/2, size(edge,2),1), edge'];
        left = [repmat(-width/2, size(edge,2),1), edge'];
        bottom = [edge' repmat(-width/2, size(edge,2),1)];
        border = cat(1, top, right, left, bottom);
    end
end

function voronoiArea = calculateVoronoiArea(uniqueKSpacePoints)
    [~, vertecies, edgeListofVolumes] = calculateVornoi(uniqueKSpacePoints);
    voronoiArea = zeros(size(edgeListofVolumes, 1), 1);
    %loop through each edge list of a voronoi hulll
    for iVolume = 1:length(edgeListofVolumes)
        %This happends when voronoi diagram goes to inf
        if(all(edgeListofVolumes{iVolume} ~= 1))
            %Calculate the volume from Voronoi diagram 
            [~, voronoiArea(iVolume)] = convhulln(vertecies(edgeListofVolumes{iVolume}, :));
        end
    end
end

function [triangulation, vertecies, edgeListofVolumes] = calculateVornoi(uniqueKSpacePoints)
    %get delaunay triangulation from points
    triangulation = delaunayTriangulation(uniqueKSpacePoints(:,1), uniqueKSpacePoints(:,2));
    
    %create voronoi diagram from triangulation
    [vertecies, edgeListofVolumes] = voronoiDiagram(triangulation);

end



function plotVoronoi(bordered_k_space, weights, indexofOccurences)
    [triangulation, kVertecies, edgeListofVolumes] = calculateVornoi(bordered_k_space);
    figure
    %plotting the sampling points
    subplot(2,2,1)
    scatter(triangulation.Points(:,1), triangulation.Points(:,2)), title('Sampling Points with Border');
    
    %plotting voronoi diagram
    subplot(2,2,2)
    hold on
    scatter(triangulation.Points(:,1), triangulation.Points(:,2), 30, '.'), title('Voronoi diagram');
    edges = cell(length(edgeListofVolumes),1);
    for i = 1:length(edgeListofVolumes)
        edges{i} = kVertecies(edgeListofVolumes{i},:);
    end
    cellfun(@(convex_hull) plot(convex_hull(:,1), convex_hull(:,2), '-r'), edges);
    hold off
    
    subplot(2,2,3)
    x = triangulation.Points(indexofOccurences, 1);
    y = triangulation.Points(indexofOccurences, 2);
    
    scatter(x, y, [], weights);
    colorbar
    
    density = 1./weights;
    subplot(2,2,4)
    x = triangulation.Points(indexofOccurences, 1);
    y = triangulation.Points(indexofOccurences, 2);
    xDiameter = max(x, [], 'all') - min(x, [], 'all');
    yDiameter = max(y, [], 'all') - min(y, [], 'all');
    xInterpVector = linspace(-xDiameter/2, xDiameter/2, length(x)/10);
    yInterpVector = linspace(-yDiameter/2, yDiameter/2, length(y)/10);
    [xMesh, yMesh] = meshgrid(xInterpVector, yInterpVector);
    zi = griddata(x, y, density, xMesh, yMesh);
    
    mesh(xMesh, yMesh, zi);
    colorbar
end

function weights = voronoi(k_space_file, NameValueArgs, k_space)
    % get border for vornoi diagrams
    border = createBorder(k_space_file, NameValueArgs.areaOfSupport);

    %get the number of unique points in the bordered_k_space. Duplicate points
    %prove problimatic for Voronoi diagrams
    [uniqueKSpacePoints, ~, indexofOccurences] = unique(k_space, 'rows', 'stable');
    
    bordered_k_space = cat(1, uniqueKSpacePoints, border);
    voronoiArea = calculateVoronoiArea(bordered_k_space);
    
    %calculate occurances of duplicate k space points
    %get occurences and unique indexes 
    [occurences, indexGroups, ~] = groupcounts(indexofOccurences);
    %rescale final_vol by the number of occurences.
    for iCoordinate = 1:length(indexGroups)
        kSpaceIndex = indexGroups(iCoordinate);
        voronoiArea(kSpaceIndex) = voronoiArea(kSpaceIndex) / occurences(iCoordinate);
    end
    
    weights = voronoiArea(indexofOccurences);
    weights = normalize(weights, 'range', [0.1, 1]);

    
    spatialPoints = calculateNumSpatialPoints(k_space_file);
    TR = length(k_space)/spatialPoints;

    weights = reshape(weights, [spatialPoints, TR]);
    if(NameValueArgs.isPlot)
        plotVoronoi(bordered_k_space, weights, indexofOccurences);
    end
end