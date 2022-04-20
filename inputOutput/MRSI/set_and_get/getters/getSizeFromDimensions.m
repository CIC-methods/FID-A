% getSizeFromDimensions.m
% Brenden Kadota, Sunnybrook 2022
%
% Description: Helper function for MRSI processing. Takes an MRSI structure and 
% cell array of dimension labels and returns their respective sizes. 
% If the dimension does not exist, or is empty, set the dimension size to 1.
%
% Input:
%   MRSIStruct = MRSI structure used in FID-A
%   dimensionLabels = cell array of strings of dimension labels (ie. {'t', 'y', 'x'})
%
% Output:
%   dimensionSizes = dimension sizes from labels. Non-existent dimensions set to 1.

function dimensionSizes = getSizeFromDimensions(MRSIStruct, dimensionLabels)
    arguments
        MRSIStruct (1, 1) struct
        dimensionLabels (1, :) cell
    end
    % get dimension numbers from labels
    dimNumbers = getMultipleDimensions(MRSIStruct, dimensionLabels);
    dataSize = getSize(MRSIStruct);

    % get zero dimensions 
    zeroIdx = dimNumbers == 0;
    biggerIndex = dimNumbers > length(dataSize);
    nonExistentIndexes = zeroIdx | biggerIndex;

    % get existing dimensions
    nonZeroIndexes = dimNumbers(~nonExistentIndexes);
    sizeOfNonZeroDimensions = dataSize(nonZeroIndexes);

    % initalize
    dimensionSizes = zeros(1, length(dimensionLabels));
    % set non existent dims to 1
    dimensionSizes(nonExistentIndexes) = 1;
    % set existing dimensions
    dimensionSizes(~nonExistentIndexes) = sizeOfNonZeroDimensions;
end
