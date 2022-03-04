function MRSIStruct = setDimsAllZero(MRSIStruct)
    dimNames = fieldnames(MRSIStruct.dims);
    for iName = 1:length(dimNames)
        MRSIStruct = setDimension(MRSIStruct, dimNames{iName}, 0);
    end
end