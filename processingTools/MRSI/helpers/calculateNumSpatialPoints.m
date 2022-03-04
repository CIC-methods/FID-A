function numSpatialPoints = calculateNumSpatialPoints(kSpaceFile)
    kTable = readtable(kSpaceFile);
    num_TR = max(kTable.TR, [], 'all');
    numSpatialPoints = height(kTable)/num_TR;
end