function newPrevPermute = removeDimPrevPermute(prevPermute, dimNumber)
    newPrevPermute = prevPermute;

    for iDim = 1:size(prevPermute, 2)
        dimensionNumber = prevPermute{1, iDim};
        if(dimensionNumber > dimNumber)
            newPrevPermute{1, iDim} = dimensionNumber - 1;
        elseif(dimensionNumber == dimNumber)
            removeIndex = iDim;
        end
    end
    newPrevPermute(:, removeIndex) = [];
end