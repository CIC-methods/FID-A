classdef  (SharedTestFixtures={matlab.unittest.fixtures.PathFixture(...
        fullfile(fileparts(mfilename('fullpath')), '..', '..', '..', 'processingTools'), ...
        'IncludeSubfolders', true)})  ...
        testop_CSIApodize < matlab.unittest.TestCase
    properties
        MRSIStruct
        kMaxX
        kMaxY
        kx
        ky

    end
    methods(TestMethodSetup)
        function createPhantom(testCase)
            MRSI = struct();
            fovX = 200;
            fovY = 200;
            fovZ = 10;
            voxelSizeX = 10;
            voxelSizeY = 10;
            voxelSizeZ = 10;

            MRSI.fov.x = fovX;
            MRSI.fov.y = fovY;
            MRSI.fov.z = fovZ;
            MRSI.voxelSize.x = voxelSizeX;
            MRSI.voxelSize.y = voxelSizeY;
            MRSI.voxelSize.z = voxelSizeZ;
            x = -fovX/2 + voxelSizeX/2:voxelSizeX:fovX/2 - voxelSizeX/2;
            y = -fovY/2 + voxelSizeY/2:voxelSizeY:fovY/2 - voxelSizeY/2;
            z = -fovZ/2 + voxelSizeZ/2:voxelSizeZ:fovZ/2 - voxelSizeZ/2;
            MRSI.coordinates.x = x;
            MRSI.coordinates.y = y;
            MRSI.coordinates.z = z;

            kxFov = 1/(x(2) - x(1));
            kyFov = 1/(y(2) - y(1));
            deltaKx = 1/fovX;
            deltaKy = 1/fovY;
            testCase.kMaxX = kxFov/2;
            testCase.kMaxY = kyFov/2;
            testCase.kx = -kxFov/2 + deltaKx/2:deltaKx:kxFov/2 - deltaKx/2;
            testCase.ky = -kyFov/2 + deltaKy/2:deltaKy:kyFov/2 - deltaKy/2;
            MRSI.data = ones(1024, 20, 20);
            MRSI.dims.t = 1;
            MRSI.flags.apodized = 0;
            testCase.MRSIStruct = MRSI;

        end
    end

    methods(Test)

        function testApodizeCosine(testCase)
            MRSILocal = testCase.MRSIStruct;
            MRSILocal.flags.spatialFT = 0;
            MRSILocal.dims.kx = 2;
            MRSILocal.dims.ky = 3;
            apodized = op_CSIApodize(MRSILocal, "functionType", 'cosine');

            W = getWeights('cosine', testCase.kx, testCase.ky, testCase.kMaxX, testCase.kMaxY);
            realData = applyApodization(MRSILocal, W, 0);
            testCase.verifyEqual(apodized.data, realData);
        end

        function testApodizeHamming(testCase)
            MRSILocal = testCase.MRSIStruct;
            MRSILocal.flags.spatialFT = 0;
            MRSILocal.dims.kx = 2;
            MRSILocal.dims.ky = 3;
            apodized = op_CSIApodize(MRSILocal, "functionType", 'hamming');

            W = getWeights('hamming', testCase.kx, testCase.ky, testCase.kMaxX, testCase.kMaxY);
            realData = applyApodization(MRSILocal, W, 0);
            testCase.verifyEqual(apodized.data, realData);

        end

        function testApodizeGaussian(testCase)
            MRSILocal = testCase.MRSIStruct;
            MRSILocal.flags.spatialFT = 0;
            MRSILocal.dims.kx = 2;
            MRSILocal.dims.ky = 3;
            apodized = op_CSIApodize(MRSILocal, "functionType", 'gaussian');

            W = getWeights('gaussian', testCase.kx, testCase.ky, testCase.kMaxX, testCase.kMaxY);
            realData = applyApodization(MRSILocal, W, 0);

            testCase.verifyEqual(apodized.data, realData);
        end

        function testApodizeCosineSpatial(testCase)
            MRSILocal = testCase.MRSIStruct;
            MRSILocal.flags.spatialFT = 1;
            MRSILocal.dims.x = 2;
            MRSILocal.dims.y = 3;
            apodized = op_CSIApodize(MRSILocal, "functionType", 'cosine');

            W = getWeights('cosine', testCase.kx, testCase.ky, testCase.kMaxX, testCase.kMaxY);
            W = fourierTransformWeights(W);
            realData = applyApodization(MRSILocal, W, 1);

            testCase.verifyEqual(apodized.data, realData, 'RelTol', 1e-10);
        end

        function testApodizeHammingSpatial(testCase)
            MRSILocal = testCase.MRSIStruct;
            MRSILocal.flags.spatialFT = 1;
            MRSILocal.dims.x = 2;
            MRSILocal.dims.y = 3;
            apodized = op_CSIApodize(MRSILocal, "functionType", 'hamming');

            W = getWeights('hamming', testCase.kx, testCase.ky, testCase.kMaxX, testCase.kMaxY);
            W = fourierTransformWeights(W);
            realData = applyApodization(MRSILocal, W, 1);

            testCase.verifyEqual(apodized.data, realData, 'RelTol', 1e-10);
        end

        function testApodizeGaussianSpatial(testCase)
            MRSILocal = testCase.MRSIStruct;
            MRSILocal.flags.spatialFT = 1;
            MRSILocal.dims.x = 2;
            MRSILocal.dims.y = 3;
            apodized = op_CSIApodize(MRSILocal, "functionType", 'gaussian');

            W = getWeights('gaussian', testCase.kx, testCase.ky, testCase.kMaxX, testCase.kMaxY);
            W = fourierTransformWeights(W);
            realData = applyApodization(MRSILocal, W, 1);

            testCase.verifyEqual(apodized.data, realData, 'RelTol', 1e-10);
        end

        function testImageAndKSpace(testCase)
            MRSILocalK = testCase.MRSIStruct;
            MRSILocalK.flags.spatialFT = 0;
            MRSILocalK.dims.kx = 2;
            MRSILocalK.dims.ky = 3;
            apodizedK = op_CSIApodize(MRSILocalK, 'functionType', 'hamming');
            apodizedK.data = fourierTransformData(apodizedK.data, apodizedK.dims.kx, apodizedK.dims.ky);
            

            MRSILocalImage = testCase.MRSIStruct;
            MRSILocalImage.flags.spatialFT = 1;

            MRSILocalImage.dims.x = 2;
            MRSILocalImage.dims.y = 3;
            MRSILocalImage.data = fourierTransformData(MRSILocalImage.data, MRSILocalImage.dims.x, MRSILocalImage.dims.y);

            apodizedImage = op_CSIApodize(MRSILocalImage, 'functionType', 'hamming');     

            testCase.verifyEqual(apodizedImage.data, apodizedK.data, 'AbsTol', 1e-10);
        end
    end
end

function realData = applyApodization(MRSILocal, W, isSpatial)
    realData = MRSILocal.data;
    realData = permute(realData, [2, 3, 1]);
    if(isSpatial)
        for iIndex = 1:size(realData, 3)
            realData(:, :, iIndex) = conv2(realData(:, :, iIndex), W, 'same');
        end
    else
        realData = realData .* W;
    end
    realData = permute(realData, [3, 1, 2]);
end


function wFT = fourierTransformWeights(W)
    if(mod(size(W, 1), 2) == 1)
        W = circshift(W, 1, 1);
    end
    if(mod(size(W, 2), 2) == 1)
        W = circshift(W, 1, 2);
    end

    wFT = fftshift(fftshift(fft2(fftshift(fftshift(W, 1), 2)), 1), 2);
    wFT = wFT/numel(wFT);
end

function W = getWeights(apodizationType, kX, kY, kMaxX, kMaxY)
    switch(apodizationType)
        case("cosine")
            W1 = cos(pi*(kX)/(2*kMaxX));
            W2 = cos(pi*(kY)/(2*kMaxY));
        case("gaussian")
            W1 = exp(-4*(kX/kMaxX).^2);
            W2 = exp(-4*(kY/kMaxY).^2);
        case('hamming')
            W1 = 0.54 + 0.46*cos(pi*kX/kMaxY);
            W2 = 0.54 + 0.46*cos(pi*kX/kMaxY);
        otherwise
            error('No function found!');
    end
    W = W1' * W2;
end

function data = fourierTransformData(data, kxDim, kyDim)
    if(mod(size(data, kxDim), 2) == 1)
        data = circshift(data, 1, kxDim);
    end
    if(mod(size(data, kyDim), 2) == 1)
        data = circshift(data, 1, kyDim);
    end

    data = fftshift(fft(fftshift(data, kxDim), [], kxDim), kxDim);
    data = fftshift(fft(fftshift(data, kyDim), [], kyDim), kyDim);
end
