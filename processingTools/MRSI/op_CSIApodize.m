% op_CSIApodize.m
%
% USAGE:
% MRSIStruct = op_CSIApodize(MRSIStruct, functionType='hamming', fwhm=[]); 
%
% DESCRIPTION:
% Applies apodization/spatial smoothing to the MRSI data. 3 spatial smoothing 
% fucntions can be used, hamming, cosine, and gaussian. Full width half maxium 
% can be passed to determine the amount of spatial smoothing.
%
% INPUTS:
% MRSIStruct         = MRSI strucuture used in FID-A
% filterArguments (name value arguments)
%   functionType     = describes filter function. Can be 'hamming', 'cosine', 'gaussian'
%   fullWidhtHalfMax = full width half maximum of filter function in mm. Default is [].
%                      This can only be used with the gaussian function.
%
% OUTPUT:
% 
function MRSIStruct = op_CSIApodize(MRSIStruct, filterArguments)
    arguments
        MRSIStruct (1,1) struct
        filterArguments.functionType (1, :) char ...
            {mustBeMember(filterArguments.functionType,...
                         {'hamming', 'cosine', 'gaussian'})} = 'hamming'
        filterArguments.fullWidthHalfMax = [];
    end
    
    kGridX = getCoordinates(MRSIStruct, 'kx');
    kGridY = getCoordinates(MRSIStruct, 'ky');
    deltaKx = kGridX(2) - kGridX(1);
    deltaKy = kGridY(2) - kGridY(1);
    kMaxX = max(abs(kGridX)) + (deltaKx)/2;
    kMaxY = max(abs(kGridY)) + (deltaKy)/2;
    switch(lower(filterArguments.functionType))
        case("cosine")
            xWeights = cos(pi*(kGridX)/(2*kMaxX));
            yWeights = cos(pi*(kGridY)/(2*kMaxY));
        case("gaussian")
            if(isempty(filterArguments.fullWidthHalfMax))
                xWeights = gaussian(kGridX, kMaxX/2);
                yWeights = gaussian(kGridY, kMaxY/2);
            else
                halfFullWidth = filterArguments.fullWidthHalfMax/2;
                sigma = sqrt(halfFullWidth^2/2*log(0.5));
                xCoordinate = getCoordinates(MRSIStruct, 'x');
                yCoordinate = getCoordinates(MRSIStruct, 'y');

                %shift back into k-space
                xWeights = gaussian(xCoordinate, sigma);
                xWeights = fftshift(fft(fftshift(xWeights)));
                yWeights = gaussian(yCoordinate, sigma);
                yWeights = fftshift(fft(fftshift(yWeights)));
            end
        case('hamming')
            xWeights = 0.54 + 0.46*cos(pi*kGridX/kMaxY);
            yWeights = 0.54 + 0.46*cos(pi*kGridX/kMaxY);
        otherwise
            error('No function found!')
    end
    weightMatrix = xWeights' * yWeights;

    
    if(getFlags(MRSIStruct, 'spatialFT') == 0)
        [MRSIStruct, prevPermute, prevSize] = reshapeDimensions(MRSIStruct, {'ky', 'kx'});
        data = getData(MRSIStruct);
        data = data.*weightMatrix;
    else
        [MRSIStruct, prevPermute, prevSize] = reshapeDimensions(MRSIStruct, {'y', 'x'});
        data = getData(MRSIStruct);
        if(mod(size(weightMatrix, 1), 2) == 1)
            weightMatrix = circshift(weightMatrix, 1, 1);
        end
        if(mod(size(weightMatrix, 2), 2) == 1)
            weightMatrix = circshift(weightMatrix, 1, 2);
        end
        weightsFT = FIDAfft(fftshift(weightMatrix, 1),1,'t');
        weightsFT = FIDAfft(fftshift(weightsFT, 2),2,'t');
        weightsFT = weightsFT/(numel(weightsFT));
        
        for iExtra = 1:getSizeFromDimensions(MRSIStruct, {'extras'})
            data(:, :, iExtra) = conv2(data(:, :, iExtra), weightsFT, 'same');
        end
    end
    MRSIStruct = setData(MRSIStruct, data);
    MRSIStruct = reshapeBack(MRSIStruct, prevPermute, prevSize);
    MRSIStruct = setFlags(MRSIStruct, 'apodized', true);

end

function y = gaussian(x, sigma)
    y = -exp(x.^2/(2*sigma.^2));
end