% reshapeDimensions.m
%
% Permute dimensions based on the input of dims. Non inputted diimensions are appeneded to the end
% INPUT
% in: MRSI strucutre
% dims: dimension array
% OUTPUT
% out: output MRSI structure with permuted dims
% prev: Previous permutation to use to permute back to original
function [MRSIStruct, pastPermute, pastSize] = reshapeDimensions(MRSIStruct, permuteLabels)

    permutationNumbers = zeros(1, length(permuteLabels));
    for iDim = 1:length(permuteLabels)
        permutationNumbers(iDim) = getDimension(MRSIStruct, permuteLabels(iDim));
    end
    
    if any(permutationNumbers == 0)
        zeroIdx = permutationNumbers == 0;
        label = permuteLabels(zeroIdx);
        error('Can not reshape with zero dimension of %s', label{:});
    end

    data = getData(MRSIStruct);
    %find dimensions not in permutationNUmbers
    remainingDimensions = setdiff(1:ndims(data), permutationNumbers);
    %create full purmutation order
    permutationOrder = [permutationNumbers, remainingDimensions];
    %permute data
    data = permute(data, permutationOrder);
    
    %get past size
    pastSize = size(data);
    dims_size = num2cell(pastSize(1:length(permutationNumbers)));
    data = reshape(data, dims_size{:}, []);
    

    %get permutation order back
    pastPermute = getPastPermuteAndLabels(permutationOrder, MRSIStruct);
    MRSIStruct = setData(MRSIStruct, data);
    MRSIStruct = updateDims(permutationNumbers, MRSIStruct, permuteLabels);
end

function pastPermute = getPastPermuteAndLabels(permutationOrder, MRSIStruct)
    pastPermute = cell(2, length(permutationOrder));
    dimMap = getMapOfDimensions(MRSIStruct);
    for iDimension = 1:length(permutationOrder)
        previousDimension = permutationOrder(iDimension);
        pastPermute{1, previousDimension} = iDimension;
        pastPermute{2, previousDimension} = dimMap(previousDimension);
    end
end

function MRSIStruct = updateDims(permutationNumbers, MRSIStruct, permuteLabels)
    MRSIStruct = setDimsAllZero(MRSIStruct);
    for iDim = 1:length(permutationNumbers)
        MRSIStruct = setDimension(MRSIStruct, permuteLabels{iDim}, iDim);
    end
    %if permutationNumber is less than the number of dimensions in the data, the
    %remaining dimensions are vectorized and put into dim field extra.
    MRSIStruct = setDimension(MRSIStruct, 'extras', length(permutationNumbers) + 1);
end
