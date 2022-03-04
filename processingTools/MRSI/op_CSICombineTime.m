function MRSIStruct = op_CSICombineTime(MRSIStruct, dimensionLabel)
    newTimeLength = prod(getSizeFromDimensions(MRSIStruct, {'t', dimensionLabel}));
    [MRSIStruct, prevPermute, prevSize]= reshapeDimensions(MRSIStruct, {'t', dimensionLabel});
    data = getData(MRSIStruct);
    
    data = reshape(data, newTimeLength, []);
    MRSIStruct = setData(MRSIStruct, data);
    prevSize(2) = [];
    prevSize(1) = newTimeLength;

    prevPermute = removeDimPrevPermute(prevPermute, 2);
    MRSIStruct = reshapeBack(MRSIStruct, prevPermute, prevSize);
end