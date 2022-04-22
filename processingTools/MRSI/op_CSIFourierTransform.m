% op_CSIFOurierTransform.m
% Brenden Kadota, SunnyBrook Hosptial 2021.
%
% Does spectral and spatial fourier transform on MRSI structure. The fast
% fourier transform is done on the spectral dimension. If the spatial k
% space is cartesian, the fast fourier transform is applied.
% USAGE:
% out = op_CSIFourierTransform(out, spatialFT, spectralFT, isCartesian)
%
%
% input:    in         = Twix object of CSI data
%           spatialFT  = Boolean flag (1 or 0) to compute fourier
%                        transfrom along spatial dimension
%           spectralFT = Same as above but for the spectral dimension
%           isCartesian= boolean flag for fft or slow fourier tranform alon
%                        the spatial dimension
%
% output:   in         = Twix object with new out.specs field of the
%                         fourier transformed data
function MRSIStruct = op_CSIFourierTransform(MRSIStruct, k_file, fourierTransform)
    arguments
        MRSIStruct (1,1) struct
        k_file (1,:) char {mustBeFileorDefault} = ""
        fourierTransform.spatial (1,1) logical {mustHaveSpatial(fourierTransform.spatial, MRSIStruct)}
        fourierTransform.spectral (1,1)  logical {mustHaveSpectral(fourierTransform.spectral, MRSIStruct)}
    end
    % set default values for spatial and spectral fourier transform
    fourierTransform = setDefaultFlags(fourierTransform, MRSIStruct);

    % spatial dimension fourier transform
    if (fourierTransform.spatial == 1)
        disp('Calculating spatial dimension');
        if(k_file == "")
            %applying the fast fourier transform if k space is cartesian
            MRSIStruct = applyFastFourierTransformSpatial(MRSIStruct);
        else
            [kTrajectory, numSpatial, numSpectral] = readKSpaceFile(k_file, MRSIStruct);
            MRSIStruct = slowFourierTransfrom(MRSIStruct, kTrajectory, numSpatial, numSpectral);
            MRSIStruct = calculateSpectralValues(MRSIStruct, numSpatial, numSpectral);
        end
        
        MRSIStruct = setFlags(MRSIStruct, 'spatialFT', true);
        %waterRemoved = op_CSIRemoveWater(MRSIStruct);
    end

    % spectral fourier transform
    if (fourierTransform.spectral)
        % fast fourier transform along the time dimension
        disp('Calculating spectral dimension');
        MRSIStruct = fastFourierTransformTime(MRSIStruct);
    end

end

function sft2_Oper = sft2_Operator(InTraj, OutTraj, Ift_flag)
    %
    % sft2 Perform 2-dimensional slow Fourier transform
    %
    % This function was written by Bernhard Strasser, April 2018.
    %
    %
    % The function performs a slow Fourier transform in two spatial dimensions by using the definition of the Fourier Transform
    % FT(f)(a_j) = sum(k=-N/2;N/2-1) f(x_k)exp(-2pi i x_k a_j)
    % (For even N. For odd, change sum slightly)
    %
    %
    % sft2Operator = sft2(A,Ifft_flag)
    %
    % Input:
    % -         Ift_flag                    ...     Flag of 0 or 1 determining if an inverse Fourier transform or normal should be performed.
    % -         InputSize                   ...     Size of Input data.
    %
    % Output:
    % -         sft2Operator                ...     Fourier transformation matrix which can be applied to the data by Out = sft2Operator*In
    %
    %
    % Feel free to change/reuse/copy the function.
    % If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
    % Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
    % File dependancy: ?
    % 0. Preparations

    % Define Exponent of exp
    if(~Ift_flag)
        Expy = -2*pi*1i;
    else
        Expy = 2*pi*1i;
    end

    % Define Output Size
    NOut = size(OutTraj,1);

    % 1. FT

    sft2_Oper = zeros([size(OutTraj,1) size(InTraj,1)]);
    for j=1:NOut
        sft2_Oper(j,:) = exp(Expy*(OutTraj(j,1)*InTraj(:,1) ...
            + OutTraj(j,2)*InTraj(:,2)));
    end

    if(Ift_flag)
        sft2_Oper = sft2_Oper / (size(InTraj,1));
    end

end

function mustHaveSpatial(a, in)
    if isfield(in, 'spatialFT') && (a == true && in.spatialFT == 1)
        eidType = 'spatialTransform:noSpatialTransform';
        msgType = 'Spaitial Fourier Transform already done!';
        throwAsCaller(MException(eidType,msgType))
    end
end

function mustHaveSpectral(a, in)
    if isfield(in, 'spectral') && (a == true && in.spectral == 1)
        eidType = 'spectralTransform:noSpectralTransform';
        msgType = 'Spectral Fourier Transform already done!';
        throwAsCaller(MException(eidType,msgType))
    end
end

function mustBeFileorDefault(file)
    if ~isfile(file) && ~strcmp(file, "")
        eidType = 'mustBeFile:notAFile';
        msgType = 'Invalid value fo k_file, must be a file';
        throwAsCaller(MException(eidType,msgType))
    end
end


function fourierTransform = setDefaultFlags(fourierTransform, in)
    if(~isfield(fourierTransform, 'spatial'))
        if(in.flags.spatialFT)
            fourierTransform.spatial = 0;
        else
            fourierTransform.spatial = 1;
        end
    end
    if(~isfield(fourierTransform, 'spectral'))
        if(in.flags.spectralFT)
            fourierTransform.spectral = 0;
        else
            fourierTransform.spectral = 1;
        end
    end
end



function MRSIStruct = applyFastFourierTransformSpatial(MRSIStruct)
    disp('Applying fast fourier transform');
    %apply half pixel shift
    MRSIStruct = halfPixelShift(MRSIStruct);
    %get data
    data = getData(MRSIStruct);
    xDimension = getDimension(MRSIStruct, 'kx');
    if(mod(getSizeFromDimensions(MRSIStruct, {'kx'}), 2) == 1)
        data = circshift(data, 1, xDimension);
    end
    data = fftshift(fft(fftshift(data, xDimension), [], xDimension), xDimension);

    yDimension = getDimension(MRSIStruct, 'ky');
    if(mod(getSizeFromDimensions(MRSIStruct,{'ky'}), 2) == 1)
        data = circshift(data, 1, yDimension);
    end
    data = fftshift(fft(fftshift(data, yDimension), [], yDimension), yDimension);
    MRSIStruct = setData(MRSIStruct, data);

    MRSIStruct = setDimension(MRSIStruct, 'x', getDimension(MRSIStruct, 'kx'));
    MRSIStruct = setDimension(MRSIStruct, 'y', getDimension(MRSIStruct, 'ky'));
    MRSIStruct = setDimension(MRSIStruct, 'z', getDimension(MRSIStruct, 'kz'));
    MRSIStruct = setDimension(MRSIStruct, 'kx', 0);
    MRSIStruct = setDimension(MRSIStruct, 'ky', 0);
    MRSIStruct = setDimension(MRSIStruct, 'kz', 0);
end


function MRSIStruct = slowFourierTransfrom(MRSIStruct, kTrajectory, numSpatial, numSpectral)
    
    [xCoordinates, yCoordinates, imageTrajectory] = getImageTrajectory(MRSIStruct);
    
    %creating fourier transform operator for spatial domain
    sftOperator = sft2_Operator(kTrajectory, imageTrajectory, 0);
    
    %permute so first 3 dimensions are x, y and t
    [MRSIStruct, prevPermute, prevSize] = reshapeDimensions(MRSIStruct, {'t', 'ky'});
    data = getData(MRSIStruct);

    % apply slow fourier transform matrix to data
    image = applySlowFourierTranformMatrix(MRSIStruct, sftOperator, data, numSpectral, numSpatial);
    MRSIStruct = setData(MRSIStruct, image);

    kyDimension = getDimension(MRSIStruct, 'ky');
    prevPermute = removeDimPrevPermute(prevPermute, kyDimension);
    prevPermute = addDimPrevPermute(prevPermute, 'y', kyDimension);
    prevPermute = addDimPrevPermute(prevPermute, 'x', kyDimension + 1);

    prevSize(1) = numSpectral;
    prevSize(2) = length(yCoordinates);
    prevSize = [prevSize(1:2), length(xCoordinates), prevSize(3:end)];
    MRSIStruct = reshapeBack(MRSIStruct, prevPermute, prevSize);
end

function image = applySlowFourierTranformMatrix(MRSIStruct, sftOperator, data, numSpectral, numSpatial)
    yLength = length(getCoordinates(MRSIStruct, 'y'));
    xLength = length(getCoordinates(MRSIStruct, 'x'));

    %image dimensions after fourier tranform is time, y, x, and extras.
    imageDimensions = [numSpectral, yLength, xLength, ...
                                        getSizeFromDimensions(MRSIStruct, {'extras'})];
    image = zeros(imageDimensions);
    
    for iPoint = 1:numSpectral
        startingPoint = (iPoint - 1)*numSpatial + 1;
        endingPoint = iPoint * numSpatial;
        kSpaceSlice = data(startingPoint:endingPoint, :, :);
        vectorizedSlice = reshape(kSpaceSlice, [], size(kSpaceSlice,3));
        ftVectorizedSlice = sftOperator*vectorizedSlice;

        imageSlice = reshape(ftVectorizedSlice, [yLength, xLength, size(ftVectorizedSlice,2)]);
        image(iPoint, :, :, :) = imageSlice;
    end
end

function [xCoordinates, yCoordinates, imageTrajectory] = getImageTrajectory(MRSIStruct)
    xCoordinates = getCoordinates(MRSIStruct, 'x');
    yCoordinates = getCoordinates(MRSIStruct, 'y');
    
    %applying the slow fourier transform if the k space is non cartesian
    [x, y] = meshgrid(xCoordinates, yCoordinates);
    imageTrajectory = [x(:), y(:)];
end

function MRSIStruct = fastFourierTransformTime(MRSIStruct)
    data = getData(MRSIStruct);
    timeDimension = getDimension(MRSIStruct, 't');
    %fourier transform in the spectral domain
    data = fftshift(fft(data, [], timeDimension), timeDimension);
    
    MRSIStruct = setData(MRSIStruct, data);
    ppm = calculatePPM(MRSIStruct);
    ppm = ppm + 4.65;
    MRSIStruct = setPPM(MRSIStruct, ppm);
    %flip ppm
    MRSIStruct = setFlags(MRSIStruct, 'spectralFT', true);
end

function ppm = calculatePPM(MRSIStruct)
    %lower bounds of frequency
    spectralWidth = getSpectralWidth(MRSIStruct);
    timeSize = getSizeFromDimensions(MRSIStruct, {'t'});

    step = spectralWidth/timeSize;
    lowerBound = -spectralWidth/2 + step/2;
    %upper bounds of frequency
    upperBound = spectralWidth/2 - step/2;
    %frequency step
    
    %calculating the frequency
    frequency=lowerBound:step:upperBound;
    
    %calculating the ppm
    ppm=-frequency/(MRSIStruct.Bo*42.577);
end


function MRSIStruct = calculateSpectralValues(MRSIStruct, numSpatial, numSpectral)
    spectralDwellTime = calculateSpectralDwellTime(MRSIStruct, numSpatial);
    spectralWidth = 1/spectralDwellTime;
    spectralTime = calculateSpectralTime(spectralDwellTime, numSpectral);

    MRSIStruct = setSpectralWidth(MRSIStruct, spectralWidth);
    MRSIStruct = setSpectralDwellTime(MRSIStruct, spectralDwellTime);
    MRSIStruct = setSpectralTime(MRSIStruct, spectralTime);
end

function spectralDwellTime = calculateSpectralDwellTime(MRSIStruct, spatialPoints)
    adcDwellTime = getAdcDwellTime(MRSIStruct);
    spectralDwellTime = spatialPoints * adcDwellTime;
end

function spectralTime = calculateSpectralTime(spectralDwellTime, spatialPoints)
    spectralTime = 0:spectralDwellTime:spectralDwellTime*(spatialPoints - 1);
end

function MRSIStruct = halfPixelShift(MRSIStruct)
    kx = getCoordinates(MRSIStruct, 'kx');
    ky = getCoordinates(MRSIStruct, 'ky');
    halfPixelX = getVoxSize(MRSIStruct, 'x')/2;
    halfPixelY = getVoxSize(MRSIStruct, 'y')/2;
    kShift = kx*halfPixelX + ky'*halfPixelY;

    [MRSIStruct, prevPermute, prevSize] = reshapeDimensions(MRSIStruct, {'ky', 'kx'});
    data = getData(MRSIStruct);
    data = data .* exp(-1i*2*pi*kShift);
    MRSIStruct = setData(MRSIStruct, data);
    MRSIStruct = reshapeBack(MRSIStruct, prevPermute, prevSize);
end
