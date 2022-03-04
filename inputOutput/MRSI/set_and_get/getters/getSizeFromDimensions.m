
function dimensionSizes = getSizeFromDimensions(MRSIStruct, dimensionLabels)
    arguments
        MRSIStruct (1, 1) struct
        dimensionLabels (1, :) cell
    end
    %get Dimension numbers
    dimNumbers = getMultipleDimensions(MRSIStruct, dimensionLabels);
    %check for zero size
    if(any(dimNumbers == 0))
        label = dimensionLabels(dimNumbers == 0);
        error('Can not get size from zero dimension %s', char(label));
    end
    dataSize = getSize(MRSIStruct);
    dataSize(end + 1:max(dimNumbers)) = 1;
    dimensionSizes = dataSize(dimNumbers);
end
