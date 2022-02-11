function MRSIStruct = removeDimension(MRSIStruct, dimLabel)
    dimNumber = getDimension(MRSIStruct, dimLabel);

    dimLabels = fieldnames(MRSIStruct.dims);
    %updating MRSI parameters
    for iName = 1:length(dimLabels)
        curDimLabel = dimLabels{iName};
        curDimNumber = getDimension(MRSIStruct, curDimLabel);
        if(curDimNumber > dimNumber)
            MRSIStruct = setDimension(MRSIStruct, curDimLabel, curDimNumber - 1);
        elseif(curDimNumber == dimNumber)
            MRSIStruct = setDimension(MRSIStruct, curDimLabel, 0);
        end
    end
end