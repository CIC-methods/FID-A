function MRSIStruct = op_CSIApodize(MRSIStruct, functions)
    arguments
        MRSIStruct (1,1) struct
        functions.functionType (1, :) char {mustBeMember(functions.functionType, {'hamming', 'cosine', 'gaussian'})} = 'hamming'
    end
    
    kGridX = getCoordinates(MRSIStruct, 'kx');
    kGridY = getCoordinates(MRSIStruct, 'ky');
    deltaKx = kGridX(2) - kGridX(1);
    deltaKy = kGridY(2) - kGridY(1);
    kMaxX = max(abs(kGridX)) + (deltaKx)/2;
    kMaxY = max(abs(kGridY)) + (deltaKy)/2;
    switch(lower(functions.functionType))
        case("cosine")
            W1 = cos(pi*(kGridX)/(2*kMaxX));
            W2 = cos(pi*(kGridY)/(2*kMaxY));
        case("gaussian")
            W1 = exp(-4*(kGridX/kMaxX).^2);
            W2 = exp(-4*(kGridY/kMaxY).^2);
        case('hamming')
            W1 = 0.54 + 0.46*cos(pi*kGridX/kMaxY);
            W2 = 0.54 + 0.46*cos(pi*kGridX/kMaxY);
        otherwise
            error('No function found!')
    end
    W = W1' * W2;

    
    if(getFlags(MRSIStruct, 'spatialFT') == 0)
        [MRSIStruct, prevPermute, prevSize] = reshapeDimensions(MRSIStruct, {'ky', 'kx'});
        data = getData(MRSIStruct);
        data = data.*W;
    else
        [MRSIStruct, prevPermute, prevSize] = reshapeDimensions(MRSIStruct, {'y', 'x'});
        data = getData(MRSIStruct);
        if(mod(size(W, 1), 2) == 1)
            W = circshift(W, 1, 1);
        end
        if(mod(size(W, 2), 2) == 1)
            W = circshift(W, 1, 2);
        end
        weightsFT = fftshift(fft(fftshift(W, 1), [], 1), 1);
        weightsFT = fftshift(fft(fftshift(weightsFT, 2), [], 2), 2);
        weightsFT = weightsFT/(numel(weightsFT));
        
        for iExtra = 1:getSizeFromDimensions(MRSIStruct, {'extras'})
            data(:, :, iExtra) = conv2(data(:, :, iExtra), weightsFT, 'same');
        end
    end
    MRSIStruct = setData(MRSIStruct, data);
    MRSIStruct = reshapeBack(MRSIStruct, prevPermute, prevSize);
    MRSIStruct = setFlags(MRSIStruct, 'apodized', true);

end
