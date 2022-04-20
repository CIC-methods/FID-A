function [kTrajectory, spatial_points, temporal_points] = readKSpaceFile(kFileName, MRSIStruct)
    kSpaceTable = readtable(kFileName);
    kTrajectory = [kSpaceTable.Kx, kSpaceTable.Ky];
    
    spatial_points = calculateNumSpatialPoints(kFileName);
    temporal_points = floor(getSizeFromDimensions(MRSIStruct, {'t'})/spatial_points);
end