classdef  (SharedTestFixtures={matlab.unittest.fixtures.PathFixture(...
        fullfile(fileparts(mfilename('fullpath')), '..', '..', '..', 'inputOutput'), ...
        'IncludeSubfolders', true)})  ...
    testio_CSIloadtwix < matlab.unittest.TestCase
    properties
        dataFile
    end
    methods(TestMethodSetup)
        function createPhantom(testCase)
            fileRoot = fileparts(mfilename('fullpath'));
            testCase.dataFile = fullfile(fileRoot, 'data', ...
                            "meas_MID00358_FID58476_csi_slaser_40_TR2000.dat");
        end
    end
    methods(Test)
        function testLoadData(testCase)
            twix = io_CSIload_twix(testCase.dataFile);
            testCase.verifyClass(twix, 'struct', 'output should be of type struct');
            
        end
        function testValues(testCase)
            twix = io_CSIload_twix(testCase.dataFile);
            testCase.verifyClass(getData(twix), 'single')
            testCase.verifyClass(getAllDims(twix), 'struct')
            testCase.verifyClass(twix.voxelSize, 'struct')
            testCase.verifyEqual(size(twix.data), getSize(twix), ['sz field should be '...
                'same size as twix.data size']);
            testCase.verifyEqual(length(twix.adcTime), getSizeFromDimensions(twix, {'t'}));
            testCase.verifyEqual(twix.adcDwellTime, twix.adcTime(2) - twix.adcTime(1));
        end
        function testError(testCase)
            testCase.verifyError(@() io_CSIload_twix(''), 'MATLAB:validators:mustBeNonzeroLengthText')
            testCase.verifyError(@() io_CSIload_twix('.'), 'MATLAB:validators:mustBeFile')
        end
    end
end

