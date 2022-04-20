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

function fig = op_CSIPlot(MRSIStruct, plotType, indecies)
    arguments
        MRSIStruct (1, 1) struct
        plotType.plane_type char {mustBeMember(plotType.plane_type, {'real', 'imag', 'abs', 'phase'})} = 'real'
        plotType.tDim char {mustBeMember(plotType.tDim, {'ppm', 't'})}
        indecies.ppmRange (1, 2) double
        indecies.coilIndex (1, 1) double = 1
        indecies.averageIndex (1, 1) double = 1
        indecies.extraIndex (1, 1) double = 1
        indecies.xIndecies (2, 1) double
        indecies.yIndecies (2, 1) double
        indecies.xRange (2, 1) double
        indecies.yRange (2, 1) double
    end
    plotType = setDefaultTimeDim(MRSIStruct, plotType);
    %check arguments to make sure they are okay
    checkArguments(MRSIStruct, indecies, plotType);
    %get x and y range from indecies
    [xRange, yRange] = getXandYRanges(indecies, MRSIStruct);

    %set default time parameters
    indecies.ppmRange = setDefaultTime(MRSIStruct, plotType.tDim, indecies);
    [ppmRange, timeVector] = getTimeBounds(MRSIStruct, indecies.ppmRange, plotType.tDim);
    
    [nonZeroIndex, data] = permuteData(MRSIStruct);
    data = indexData(ppmRange, xRange, yRange, indecies, nonZeroIndex, data);
    data = changePlane(plotType, data);
    

    [scalefactorX, scalefactorY] = getScaleFactors(data, MRSIStruct, timeVector);

    xCoordinate = getCoordinates(MRSIStruct, 'x');
    voxSizeX = getVoxSize(MRSIStruct, 'x');
    yCoordinate = getCoordinates(MRSIStruct, 'y');
    voxSizeY = getVoxSize(MRSIStruct, 'y');
    %start time vector at zero
    if(getFlags(MRSIStruct, 'spectralFT') == 1)
         %timeVector = flip(timeVector);
         timeVector = flip(timeVector);
    end
    timeVector = timeVector - timeVector(1);
    %scale to be plotted in a voxel 
    timeVector = timeVector * scalefactorX;
    %offset ppm by starting coordinate
    timeVector = timeVector - voxSizeX/2 + xCoordinate(xRange(1));

    data = data.*scalefactorY;

    %create figure and hold
    fig = figure;
    ax = axes;
    hold(ax, 'on');

    xlim(ax, [xCoordinate(xRange(1)) - voxSizeX, xCoordinate(xRange(2)) + voxSizeX]);
    xticks(ax, xCoordinate(xRange(1):xRange(2)));
    ylim(ax, [yCoordinate(yRange(1)) - voxSizeY, yCoordinate(yRange(2)) + voxSizeY]);
    yticks(ax, yCoordinate(yRange(1):yRange(2)));

    for x_idx = 1:size(data, 3)
        for y_idx = 1:size(data, 2)
            %first scale the ppm scale so that the range is correct;
            timeVectorPlot = timeVector + (x_idx - 1)*voxSizeX;
            y_coords = data(:, y_idx, x_idx) + (size(data, 2) - y_idx)*voxSizeY + yCoordinate(yRange(1));
            %Now start plotting
            plot(ax, timeVectorPlot, y_coords, 'PickableParts', 'none', 'LineWidth', 1);
        end
    end
    hold(ax, 'off');
end

function checkArguments(in, indecies, plotType)
    %Argument checks
    coilDimension = getDimension(in, 'coils');
    if(coilDimension == 0 && indecies.coilIndex > 1)
        error('op_CSIPlot:argumentError', 'No coils to plot')
    end
    if(coilDimension ~= 0 && in.sz(coilDimension) < indecies.coilIndex)
        error('op_CSIPlot:argumentError', 'Coil index is larger than number of coils')
    end
    
    averageDimension = getDimension(in, 'averages');
    if(averageDimension == 0 && indecies.averageIndex > 1)
        error('op_CSIPlot:argumentError', 'No averages to plot')
    end
    if(averageDimension ~= 0 && in.sz(averageDimension) < indecies.averageIndex)
        error('op_CSIPlot:argumentError', 'Average index is larger than number of averages')
    end

    if(exist('indecies.x_indecies', 'var') && exist('indecies.x_range', 'var'))
        error('op_CSIPlot:argumentError', 'Only x_indecies or x_range should be used, not both');
    end
    if(exist('indecies.y_indecies', 'var') && exist('indecies.y_range', 'var'))
        error('op_CSIPlot:argumentError', 'Only y_indecies or y_range should be used, not both');
    end
    if(strcmp(plotType.tDim, 'ppm') && ~isfield(in, 'ppm'))
        error('Plotting spectrum but fids exist. Fourier transform along the spectral dimension');
    end
end

function [x, y] = getXandYRanges(indecies, in)
    xCoordinates = getCoordinates(in, 'x');
    if(~isfield(indecies, 'xIndecies') && ~isfield(indecies, 'xRange'))
        x = [1, length(xCoordinates)];
    elseif(isfield(indecies, 'xRange'))
        x_val = xCoordinates > indecies.xRange(1) & xCoordinates < indecies.xRange(2);
        x = [find(x_val, 1), find(x_val, 1, 'last')];
    else
        x = indecies.xIndecies;
    end


    yCoordinates = getCoordinates(in, 'y');
    if(~isfield(indecies, 'yIndecies') && ~isfield(indecies, 'yRange'))
        y = [1, length(yCoordinates)];
    elseif(isfield(indecies, 'yRange'))
        y_val = yCoordinates > indecies.yRange(1) & yCoordinates < indecies.yRange(2);
        y = [find(y_val, 1), find(y_val, 1, 'last')];
    else
        y = indecies.yIndecies;
    end
end

function [timeBounds, timeVector] = getTimeBounds(MRSIStruct, timeArguments, timeType)
    minTime = timeArguments(1);
    maxTime = timeArguments(2);
    if(strcmp(timeType, 't'))
        time = getSpectralTime(MRSIStruct);
        timeBounds = minTime <= time & maxTime >= time;
        timeVector = time(timeBounds);
    elseif(strcmp(timeType, 'ppm'))
        ppm = getPPM(MRSIStruct);
        timeBounds = minTime <= ppm & maxTime >= ppm;
        timeVector = ppm(timeBounds);
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

function data = indexData(timeRange, xRange, yRange, indecies, nonZeroIndex, data)
    %index for the data
    index = {timeRange, yRange(1):yRange(2), xRange(1):xRange(2), ...
        indecies.coilIndex, indecies.averageIndex, indecies.extraIndex};
    %remove dimensions that are not exsitent
    
    index = index(nonZeroIndex);
    %index data
    data = data(index{:});
end


function data = changePlane(plotType, data)
    %get lowercase
    plottingPlane = lower(plotType.plane_type);
    if(strcmp(plottingPlane, 'real'))
        data = real(data);
    elseif (strcmp(plottingPlane, 'imag'))
        data = imag(data);
    elseif (strcmp(plottingPlane, 'abs'))
        data = abs(data);
    else
        dataSize = size(data);
        for i = 1:prod(dataSize(2:end))
            data(:, i) = phase(data(:, i));
        end
    end
end

function timeRange = setDefaultTime(MRSIStruct, plotType, indecies)
    if(isfield(indecies, 'ppmRange'))
        timeRange = indecies.ppmRange;
    else
        if strcmp(plotType, 't')
            time = getSpectralTime(MRSIStruct);
            minTime = min(time);
            maxTime = max(time);
        elseif(strcmp(plotType, 'ppm'))
            ppm = getPPM(MRSIStruct);
            minTime = min(ppm);
            maxTime = max(ppm);
        end
        timeRange = [minTime, maxTime];
    end
end

function [scalefactorX, scalefactorY] = getScaleFactors(data, MRSIStruct, timeVector)
    %max spectrum difference in a voxel
    spectrumHight = max(max(data, [], 1) - min(data, [], 1), [], 'all');

    %scale factors to fit at each (x,y) coordinates
    
    voxSizeX = getVoxSize(MRSIStruct, 'x');
    voxSizeY = getVoxSize(MRSIStruct, 'y');
    scalefactorX = abs((0.8 * voxSizeX) / (timeVector(end) - timeVector(1)));
    scalefactorY = abs((0.8 * voxSizeY) / spectrumHight);
end

function plotType = setDefaultTimeDim(MRSIStruct, plotType)
    if(~isfield(plotType, 'tDim'))
        if(getFlags(MRSIStruct, 'spectralft'))
            plotType.tDim = 'ppm';
        else
            plotType.tDim = 't';
        end
    end
end