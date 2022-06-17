function lipidMask = op_CSICreateLipidMask(MRSIStruct)
    checkArguments(MRSIStruct);
    dimensionIndexes = checkDimensions(MRSIStruct);

    
end

function checkArguments(MRSIStruct) 
    checkSpatialFT(MRSIStuct);
    checkSpectralFT(MRSIStuct);
end

function [dimensionIndexs] = checkDimensions(MRSIStruct) 
    dataSize = getSize(MRSIStruct);
    dims = getDimension(MRSIStruct);
    
end
     
