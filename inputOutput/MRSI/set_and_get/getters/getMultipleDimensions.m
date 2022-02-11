function dimNumber = getMultipleDimensions(MRSIStruct, dimLabelCell)
    arguments
        MRSIStruct (1,1) struct
        dimLabelCell (1,:) cell
    end
    
    dimNumber = zeros(1, length(dimLabelCell));
    for iLabel = 1:length(dimLabelCell)
        dimNumber(iLabel) = getDimension(MRSIStruct, dimLabelCell{iLabel});
    end
end
