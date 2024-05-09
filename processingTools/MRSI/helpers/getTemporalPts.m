% getTemporalPts.m
% Jamie Near, Sunnybrook Research Institute 2024.
%
% USAGE: 
% NPtemporal = getTemporalPts(kTable,MRSIStruct)
%
% DESCRIPTION:
%This function accepts as input a "kTable" (the output of the function 
%'readkFile.m' for a rapid MRSI seuqence and the corresponding MRSI 
%structure 'MRSIStruct' and returns the number of acquired temporal points 
%per k-space position.  This is equivalent to the number of points in the 
%FID, or the number of spectral points.   
%
%INPUTS:
% kTable                 = Matlab table specifying the k-space coordinates 
%                          of a rapid MRSI acquisition.  This kTable can be
%                          obtained by reading a k-space coordinates text
%                          file using the function 'readKFile.m'.  
% MRSIStruct             = MRSI data structure in FID-A format (prior to
%                          FT). 
%
%OUTPUTS
% NPtemporal:   = The calculated number of acquired temporal points per
%                 k-space position.  


function NPtemporal = getTemporalPts(kTable,MRSIStruct)

    kPtsPerCycle = getKPtsPerCycle(kTable);

    %Find the number of temporal points:
    NPtemporal = floor(getSizeFromDimensions(MRSIStruct,{'t'})/kPtsPerCycle);

end