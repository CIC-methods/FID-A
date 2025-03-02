% op_CSIplot.m
% Brenden Kadota, McGill University 2019.
%
% USAGE:
% sim_plotCSI(in)
% sim_plotCSI(in, 'min_range', 0, 'max_range', 4, 'coilNum', 3, 'averageIndex', 1,
%             'x_indecies', [10,20], 'y_range', [-100, 100])
%
% DESCRIPTION:
% This function takes a processed MRSI twix file and plots the data using
% the matlab plot function. Horizontal coordinates are the x coordinates
% and vertical are y coordinates. The plots at specific (x,y) positions are
% the spectral coordintaes
%
% INPUTS:
% in          = input cell array of simulated spectra from a spatially resolved simulation
%   NAME VALUE ARGUMENTS (ALL OPTIONAL)
%   Name: 'min_range' Value: double, minimum ppm range or time range to plot
%   Name: 'max_range' Value: double, maximum ppm range or time range to plot
%   Name: 'coilNum' Value: integer, index of coil to plot
%   Name: 'averageIndex' Value: integer, index of average to plot
%   Name: 'x_indecies' Value: (2,1) integer, bounds (inclusive) of x indecies to plot.
%                              First index must be smaller than the second.
%   Name: 'y_indecies' Value: (2,1) integer, bounds (inclusive) of y indecies to plot.
%                              First index must be smaller than the second.
%   Name: 'x_range' Value: (2,1) double, bounds (inclusive) of x values to plot.
%                              First index must be smaller than the second.
%   Name: 'y_range' Value: (2,1) double, bounds (inclusive) of y indecies to plot.
%                              First index must be smaller than the second.
%
%
%
% OUTPUTS:
% fig   = figure handle.

function [fig] = op_CSIPlot(MRSIStruct, plotType, dimensionIndexes)
    arguments
        MRSIStruct (1, 1) struct
        plotType.plane_type char {mustBeMember(plotType.plane_type, {'real', 'imag', 'abs', 'phase'})} = 'real'
        dimensionIndexes.ppmBounds (1, 2) double
        dimensionIndexes.coilIndex (1, 1) double = 1
        dimensionIndexes.averageIndex (1, 1) double = 1
        dimensionIndexes.extraIndex (1, 1) double = 1
        dimensionIndexes.xIndecies (2, 1) double
        dimensionIndexes.yIndecies (2, 1) double
        dimensionIndexes.xRange (2, 1) double
        dimensionIndexes.yRange (2, 1) double
        dimensionIndexes.yMul (1, 1) double = 1
        dimensionIndexes.lineWidth (1, 1) double = 1
    end
    % check arguments to make sure they are okay
    checkArguments(MRSIStruct, dimensionIndexes);

    % get x and y range from indecies
    [xRange, yRange] = setDefaultXandYBounds(MRSIStruct, dimensionIndexes);
    [ppmRange, ppmVector] = getTimeBounds(MRSIStruct, dimensionIndexes);
   
    data = indexData(ppmRange, xRange, yRange, dimensionIndexes, MRSIStruct);
    data = changePlaneToRealImagOrAbs(plotType, data);
    ppmVector = setPPMVectorStartToZero(MRSIStruct, ppmVector);
    %scale to be plotted in a voxel 
    %offset ppm by starting coordinate
    [ppmVector, data] = scalePPMandData(ppmVector, data, MRSIStruct, dimensionIndexes);

    xCoordinate = getCoordinates(MRSIStruct, 'x');
    voxSizeX = getVoxSize(MRSIStruct, 'x');
    yCoordinate = getCoordinates(MRSIStruct, 'y');
    voxSizeY = getVoxSize(MRSIStruct, 'y');
    %create figure and hold
    [fig, ax] = setUpFigureAndAxes(xRange, yRange, MRSIStruct);
    
    for x_idx = 1:size(data, 3)
        for y_idx = 1:size(data, 2)
            %first scale the ppm scale so that the range is correct;
            timeVectorPlot = ppmVector + (x_idx - 1)*voxSizeX - voxSizeX/2 + xCoordinate(xRange(1));
            % scale y vector to the right position
            y_coords = data(:, y_idx, x_idx) + (size(data, 2) - y_idx)*voxSizeY + yCoordinate(yRange(1));
            %Now start plotting
            plot(ax, timeVectorPlot, y_coords, 'PickableParts', 'none', 'LineWidth', dimensionIndexes.lineWidth);
        end
    end
    hold(ax, 'off');
end



function [xBounds, yBounds] = setDefaultXandYBounds(MRSIStruct, indecies)
    xBounds = getSpatialBoundsFromIndeciesOrCoords(MRSIStruct, indecies, 'x');
    yBounds = getSpatialBoundsFromIndeciesOrCoords(MRSIStruct, indecies, 'y');
end

function [ppmLogicalIndex, ppmVector] = getTimeBounds(MRSIStruct, indecies)
    % set default ppm range
    [ppmBounds, isPPMDimension] = getDefaultPPMBounds(MRSIStruct, indecies);
    
    minPPM = ppmBounds(1);
    maxPPM = ppmBounds(2);
    
    if(isPPMDimension)
        ppm = getPPM(MRSIStruct);
        ppmLogicalIndex = ppm >= minPPM & ppm <= maxPPM;
        ppmVector = ppm(ppmLogicalIndex);
    else
        time = getSpectralTime(MRSIStruct);
        timeLogicalIndex = time >= minPPM & time <= maxPPM;
        timeVector = time(timeLogicalIndex);
        ppmLogicalIndex = timeLogicalIndex;
        ppmVector = timeVector;
    end
end

function [indexs, data] = permuteData(MRSIStruct)
    %permute data into a resonable order
    dimensionOrder = {'t', 'y', 'x', 'coils', 'averages', 'extras'};
    dimensionOrderKSpace = {'t', 'ky', 'kx', 'coils', 'averages', 'extras'};
    dimNumbers = getMultipleDimensions(MRSIStruct, dimensionOrder);
    dimNumbersKSpace = getMultipleDimensions(MRSIStruct, dimensionOrderKSpace);
    finalDim = zeros(size(dimNumbers));
    for iDim = 1:length(dimNumbers)
        if(dimNumbers(iDim) ~= 0)
            finalDim(iDim) = dimNumbers(iDim);
        end
        if(dimNumbersKSpace(iDim) ~= 0)
            finalDim(iDim) = dimNumbersKSpace(iDim);
        end
    end
    indexs = dimNumbers > 0 | dimNumbersKSpace > 0;
    data = getData(MRSIStruct);
    data = permute(data, finalDim(indexs));
end

function data = indexData(timeRange, xRange, yRange, indecies, MRSIStruct)
    [nonZeroIndex, data] = permuteData(MRSIStruct);
    %index for the dataindexData
    index = {timeRange, yRange(1):yRange(2), xRange(1):xRange(2), ...
        indecies.coilIndex, indecies.averageIndex, indecies.extraIndex};
    %remove dimensions that are not exsitent
    index = index(nonZeroIndex);
    %index data
    data = data(index{:});
end


function data = changePlaneToRealImagOrAbs(plotType, data)
    %get lowercase
    plottingPlane = lower(plotType.plane_type);
    if(strcmp(plottingPlane, 'real'))
        data = real(data);
    elseif (strcmp(plottingPlane, 'imag'))
        data = imag(data);
    elseif (strcmp(plottingPlane, 'abs'))
        data = abs(data);
    end
end

% set default time range 
function [timeBounds, isPPMDimension] = getDefaultPPMBounds(MRSIStruct, indecies)
    % if indecies has ppmRange value don't set default
    if(isfield(indecies, 'ppmBounds'))
        isPPMDimension = true;
        timeBounds = indecies.ppmBounds;
        return
    end   
    if(getFlags(MRSIStruct, 'spatialft'))
        isPPMDimension = true;
    else
        isPPMDimension = false;
    end
    if isPPMDimension
        ppm = getPPM(MRSIStruct);
    else
        ppm = getSpectralTime(MRSIStruct);
    end
    minPPM = min(ppm);
    maxPPM = max(ppm);
    % set time range indecies
    timeBounds = [minPPM, maxPPM];
end

function [scalefactorX, scalefactorY] = getPlottingScaleFactors(data, MRSIStruct, timeVector, dimensionIndexes)

    if nargin<4
        yMul=1;
    else
        yMul=dimensionIndexes.yMul;
    end

    %max spectrum difference in a voxel
    spectrumHight = max(max(data, [], 1) - min(data, [], 1), [], 'all');

    %scale factors to fit at each (x,y) coordinates
    
    voxSizeX = getVoxSize(MRSIStruct, 'x');
    voxSizeY = getVoxSize(MRSIStruct, 'y');
    scalefactorX = abs((0.8 * voxSizeX) / (timeVector(end) - timeVector(1)));
    scalefactorY = yMul * abs((0.8 * voxSizeY) / spectrumHight);
end

% Argument checks
function checkArguments(in, indecies)
    % check indexing
    checkDimensionIndexingforErrors(in, indecies.coilIndex, 'coils');
    checkDimensionIndexingforErrors(in, indecies.averageIndex, 'averages');
    
    % only indecies or rage can be accepted
    if(exist('indecies.x_indecies', 'var') && exist('indecies.x_range', 'var'))
        error('op_CSIPlot:argumentError', 'Only x_indecies or x_range should be used, not both');
    end
    if(exist('indecies.y_indecies', 'var') && exist('indecies.y_range', 'var'))
        error('op_CSIPlot:argumentError', 'Only y_indecies or y_range should be used, not both');
    end
end

function checkDimensionIndexingforErrors(in, dimensionIndex, dimensionLabel)
    % check if dimension exists and dimension is being indexed
    dimensionNumber = getDimension(in, dimensionLabel);
    if(dimensionNumber == 0 && dimensionIndex > 1)
        error('op_CSIPlot:argumentError', '%s index is larger than 1! No %s dimension exist', ...
                                            dimensionLabel, dimensionLabel)
    end
    % check if plotting index is larger than dimension index
    dimensionSize = getSizeFromDimensions(in, {dimensionLabel});
    if(dimensionIndex > dimensionSize)
        error('op_CSIPlot:argumentError', '%s index is larger than number of coils', dimensionLabel)
    end
end

% take coordinate index field and get index bounds (ie. [2, 6] is second till sixth index to plot). 
function coordinateBounds = getSpatialBoundsFromIndeciesOrCoords(in, indecies, dimensionLabel)
    Coordinates = getCoordinates(in, dimensionLabel);
    
    indeciesField = [dimensionLabel, 'Indecies'];
    rangeField = [dimensionLabel, 'Range'];

    if(isfield(indecies, indeciesField))
        coordinateBounds = indecies.xIndecies;
    elseif(isfield(indecies, rangeField))
        coordinateRange = indecies.(rangeField);
        x_val = Coordinates > coordinateRange(1) & Coordinates < coordinateRange(2);
        coordinateBounds = [find(x_val, 1), find(x_val, 1, 'last')];
    else
        coordinateBounds = [1, length(Coordinates)];
    end
end

% scale down ppm and data scale so that they fit into one voxel worth of
% plotting in MATLAB. 
function [ppmVector, data] = scalePPMandData(ppmVector, data, MRSIStruct, dimensionIndexes)
    [scalefactorX, scalefactorY] = getPlottingScaleFactors(data, MRSIStruct, ppmVector, dimensionIndexes);
    ppmVector = ppmVector * scalefactorX;
    data = data .* scalefactorY;
end


% Set the scale of PPM vector to start at zero
function ppmVector = setPPMVectorStartToZero(MRSIStruct, ppmVector)
    %start time vector at zero
    if(getFlags(MRSIStruct, 'spectralFT') == 1)
        %timeVector = flip(timeVector);
        ppmVector = flip(ppmVector);
    end
    ppmVector = ppmVector - ppmVector(1);
end



% set up figure and axes limits and ticks
function [fig, ax] = setUpFigureAndAxes(xRange, yRange, MRSIStruct)
    xCoordinate = getCoordinates(MRSIStruct, 'x');
    yCoordinate = getCoordinates(MRSIStruct, 'y');
    voxelSizeX = getVoxSize(MRSIStruct, 'x');
    voxelSizeY = getVoxSize(MRSIStruct, 'y');

    fig = figure;
    ax = axes;
    hold(ax, 'on');
    
    xlim(ax, [xCoordinate(xRange(1)) - voxelSizeX, xCoordinate(xRange(2)) + voxelSizeX]);
    xticks(ax, xCoordinate(xRange(1):xRange(2)));

    ylim(ax, [yCoordinate(yRange(1)) - voxelSizeY, yCoordinate(yRange(2)) + voxelSizeY]);
    yticks(ax, yCoordinate(yRange(1):yRange(2)));
end

