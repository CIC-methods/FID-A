% op_B0Correction.m
% Jamie Near, McGill University 2021.
%E dits from Brenden Kadota, 2021.
%
% USAGE:
% [out, phaseMap, freqMap]=op_B0Correction(in, (optional) phaseMap);
%
% DESCRIPTION:
% Corrects for slight B0 shifts by aligning all voxels to the voxel with
% the highest water peak. (Possibly not correct output as highest peak may
% no have the correct ppm shifts)
%
% INPUTS:
% in   = MRSI struct
% phaseMap = phaseMap matrix used to apply corrections
%
% OUTPUTS:
% B0Map        = 2D array of shift intensity at each specific voxel
% phaseMap     = 2D array of phase shift intensity at each specific voxel
% freqMap      = 2D array of frequency shift intensity at each specific voxel

function [MRSIStruct, phaseMap, freqMap] = op_CSIB0Correction(MRSIStruct, phaseMap, ...
        freqMap, plottingArguments)
    arguments
        MRSIStruct (1, 1) struct
        phaseMap double = []
        freqMap double = []
        plottingArguments.isPlot = false;
    end
    %only correct if coils have be combined
    if(getFlags(MRSIStruct, 'addedrcvrs') == 0)
        error("MRSIStruct Error. Add coils first in MRSIStruct");
    end
    if(getFlags(MRSIStruct, 'spatialft') == 0)
        error("Please fourier transform along the spatial dimension");
    end
    if(getFlags(MRSIStruct, 'spectralft') == 0)
        error("Please fourier transform along the spectral dimension");
    end

    %find the coordinate of the highest intensity
    [x, y, a] = findMaxCoord(MRSIStruct);

    %create twix object from coordinate with highest intensity
    referenceMRS = op_CSItoMRS(MRSIStruct, x, y, "averageIndex", a);
    referenceMRS = setFlags(referenceMRS, 'averaged', true);
    referenceMRS.flags.isISIS = 0;


    if(~isempty(phaseMap) && ~isempty(freqMap))
        [MRSIStruct, prevPermute, prevSize] = reshapeDimensions(MRSIStruct, {'t', 'y', 'x'});
        MRSIStruct = applyFreqMap(MRSIStruct, freqMap);
        MRSIStruct = applyPhaseMap(MRSIStruct, phaseMap);
        MRSIStruct = reshapeBack(MRSIStruct, prevPermute, prevSize);

    else

        [MRSIStruct, prevPermute, prevSize] = reshapeDimensions(MRSIStruct, {'t', 'y', 'x'});

        phaseMap = zeros(getSizeFromDimensions(MRSIStruct, {'y', 'x', 'extras'}));
        freqMap = zeros(getSizeFromDimensions(MRSIStruct, {'y', 'x', 'extras'}));

        data = getData(MRSIStruct);
        for e = 1:getSizeFromDimensions(MRSIStruct, {'extras'})
            for y = 1:getSizeFromDimensions(MRSIStruct, {'y'})
                for x = 1:getSizeFromDimensions(MRSIStruct, {'x'})
                    %get twix object at x and y coordinate
                    alignMRS = op_CSItoMRS(MRSIStruct, x, y, 'Extra', e);
                    alignMRS = setFlags(alignMRS, 'averaged', true);
                    alignMRS.flags.isISIS = 0;
                    %aligns twix object to the maximum intensity twix object
                    [alignMRS, phase, freq] = op_alignScans(alignMRS, referenceMRS, referenceMRS.t(end));

                    timeDimension = getDimension(alignMRS, 't');
                    %adds the aligned twix to the specs of MRSI object.
                    data(:, y, x, e) = flip(alignMRS.specs, timeDimension);
                    %add phase and frequency to maps
                    phaseMap(y, x, e) = phase;
                    freqMap(y, x, e) = freq;
                end
            end
        end
        MRSIStruct = setData(MRSIStruct, data);
        MRSIStruct = reshapeBack(MRSIStruct, prevPermute, prevSize);

        %plot maps
        if(plottingArguments.isPlot)
            plotFreqAndPhaseMap(phaseMap, MRSIStruct, freqMap);
        end
        MRSIStruct = setFlags(MRSIStruct, 'phasecorrected', 1);
        MRSIStruct = setFlags(MRSIStruct, 'freqcorrected', 1);
    end
end


%finds the (x,y) coordinate with the highest intensity spectrum.
function [x,y,a] = findMaxCoord(in)
    [~, i] = max(getData(in),[], 'all', 'linear');
    out = cell(size(in.sz));
    [out{:}] = ind2sub(in.sz, i);
    x = out{in.dims.x};
    y = out{in.dims.y};
    if(in.dims.averages)
        a = out{in.dims.averages};
    else
        a = 1;
    end


end

function plotFreqAndPhaseMap(phaseMap, MRSIStruct, freqMap)
    for i = 1:size(phaseMap,3)
        figure;
        subplot(2,1,1);
        imagesc(MRSIStruct.coordinates.x, MRSIStruct.coordinates.y, phaseMap(:,:,i)')
        title("phaseMap")
        subplot(2,1,2);
        imagesc(MRSIStruct.coordinates.x, MRSIStruct.coordinates.y, freqMap(:,:,i)')
        title("freqMap")
    end
end

function MRSIStruct = applyFreqMap(MRSIStruct, freqMap)
    %applying the freqMap map if added as an argument
    for e = 1:getSizeFromDimensions(MRSIStruct, {'extras'})
        for x = 1:getSizeFromDimensions(MRSIStruct, {'x'})
            for y = 1:getSizeFromDimensions(MRSIStruct, {'y'})
                %create twix object from the coordinates
                alignMRS = op_CSItoMRS(MRSIStruct, x, y, 'extraIndex', e);
                %align twix object with op_freqshift to the equivalent
                %freqmap coordinate.
                alignMRS = op_freqshift(alignMRS, freqMap(x, y));
                %add the specs to MRSI object
                MRSIStruct.data(:, x, y, e) = alignMRS.specs;
            end
        end
    end
end

function MRSIStruct = applyPhaseMap(MRSIStruct, phaseMap)
    %applying the freqMap map if added as an argument
    for e = 1:getSizeFromDimensions(MRSIStruct, {'extras'})
        for x = 1:getSizeFromDimensions(MRSIStruct, {'x'})
            for y = 1:getSizeFromDimensions(MRSIStruct, {'y'})
                %create twix object from the coordinates
                alignMRS = op_CSItoMRS(MRSIStruct, x, y, 'extraIndex', e);
                %align twix object with op_freqshift to the equivalent
                %freqmap coordinate.
                alignMRS = op_addphase(alignMRS, phaseMap(x, y));
                %add the specs to MRSI object
                MRSIStruct.data(:, x, y, e) = alignMRS.specs;
            end
        end
    end
end
