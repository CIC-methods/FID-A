function [MRSIStruct, phaseMap] = op_CSI0thOrderPhase(MRSIStruct, referencePoint)
arguments
    MRSIStruct (1, 1) struct
    referencePoint (1, 1) double = 1
end
    checkArguments(MRSIStruct);
    [MRSIStruct, prevPermute, prevShape] = reshapeDimensions(MRSIStruct, {'t', 'y', 'x'});
    data = getData(MRSIStruct);
    phaseMap = zeros(getSizeFromDimensions(MRSIStruct, {'y', 'x', 'extras'}));
    for e = 1:getSizeFromDimensions(MRSIStruct, {'extras'})
        for x = 1:getSizeFromDimensions(MRSIStruct, {'x'})
            for y = 1:getSizeFromDimensions(MRSIStruct, {'y'})
                mrs = op_CSItoMRS(MRSIStruct, x, y, 'extraIndex', e); 
                firstOrderPhase = 180*angle(mrs.fids(referencePoint))/pi;
                mrs = op_addphase(mrs, -firstOrderPhase, 0, 0, 1);
                phaseMap(y, x, e) = firstOrderPhase;
                timeDimension = getDimension(MRSIStruct, {'t'});
                if(getFlags(MRSIStruct, 'spectralFT') == 1)
                    data(:, y, x, e) = flip(mrs.specs, timeDimension); 
                else
                    data(:, y, x, e) = mrs.fids;
                end

            end
        end
    end
    MRSIStruct = setData(MRSIStruct, data);
    MRSIStruct = reshapeBack(MRSIStruct, prevPermute, prevShape);

end

function checkArguments(MRSIStruct)
    if(getFlags(MRSIStruct, 'spatialFT') == 0)
        error('Please fourier transform along the spectral dimension')
    end
end
