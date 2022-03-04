function newPrevPermute = addDimPrevPermute(prevPermute, dimLabel, dimNumber)
    newPrevPermute = prevPermute;

    for iDim = 1:size(prevPermute, 2)
        dimensionNumber = prevPermute{1, iDim};
        if(dimensionNumber >= dimNumber)
            newPrevPermute{1, iDim} = dimensionNumber + 1;
        end
    end
    newIndex = size(prevPermute, 2) + 1;
    newPrevPermute{1, newIndex} = dimNumber;
    newPrevPermute{2, newIndex} = dimLabel;
end