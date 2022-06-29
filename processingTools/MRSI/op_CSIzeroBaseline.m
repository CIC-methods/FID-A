% op_CSIZeroBaseline
%
% Takes mean of ppm baseline and substracts spectra by this baseline. This makes
% the ppm baseline zero.
%
% INPUT
% MRSIStruct            = MRSI structure used in FID-A
% ppmBaselineRange      = ppm range of the baseline
%
% OUTPUT
% MRSIStruct            = MRSI ouput

function MRSIStruct = op_CSIzeroBaseline(MRSIStruct, ppmBaselineRange)
    [MRSIStruct, prevPermute, prevShape] = reshapeDimensions(MRSIStruct, {'t', 'y', 'x'});
    data = getData(MRSIStruct);
    ppm = getPPM(MRSIStruct);
    ppmBaselineIndex = ppm > ppmBaselineRange(1) & ppm < ppmBaselineRange(2);
    for e = 1:getSizeFromDimensions(MRSIStruct, {'extras'})
        for x = 1:getSizeFromDimensions(MRSIStruct, {'x'})
            for y = 1:getSizeFromDimensions(MRSIStruct, {'y'})
                voxelBaselineData = mean(data(ppmBaselineIndex, y, x, e));
                realBaseline = real(voxelBaselineData);
                imagBaseline = imag(voxelBaselineData);
                data(:, y, x, e) = data(:, y, x, e) - realBaseline - 1i * imagBaseline;
            end
        end
    end
    MRSIStruct = setData(MRSIStruct, data);
    MRSIStruct = reshapeBack(MRSIStruct, prevPermute, prevShape);

end