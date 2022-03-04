function dimensionMap = getMapOfDimensions(MRSIStruct)
    values = cell2mat(struct2cell(MRSIStruct.dims));
    dimsFieldNames = fieldnames(MRSIStruct.dims);
    zeroIndexs = values == 0;
    nonZeroFieldName = dimsFieldNames(~zeroIndexs);
    dimensionMap = containers.Map(values(~zeroIndexs), nonZeroFieldName);
end