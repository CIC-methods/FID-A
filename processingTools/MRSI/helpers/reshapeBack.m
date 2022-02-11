% Used to permute back from pervious permutation from permute_dims
% INPUT
% in: MRSI structure that was permuted
% prev: prev array from premute_dims
% 
% OUTPUT
% out: permuted MRSI structure
function MRSIStruct = reshapeBack(MRSIStruct, prevPermute, prevSize)
    data = getData(MRSIStruct);
    data = reshape(data, prevSize);

    permuteBackOrder = cell2mat(prevPermute(1, :));
    permuteBackLabels = prevPermute(2, :);
    data = permute(data, permuteBackOrder);

    MRSIStruct = setDimsAllZero(MRSIStruct);
    for iDim = 1:length(permuteBackLabels)
        MRSIStruct = setDimension(MRSIStruct, permuteBackLabels{iDim}, iDim);
    end
    
    MRSIStruct = setData(MRSIStruct, data);

end
