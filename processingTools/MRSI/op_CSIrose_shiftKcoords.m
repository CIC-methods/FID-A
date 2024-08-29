% op_CSIrose_shiftKcoords.m
% Jamie Near, Sunnybrook Research Institute 2024.
%
% USAGE:
% [shiftedKTable, ind]=op_CSIrose_shiftKcoords(kFile, in);
%
% DESCRIPTION:
%This function adjusts the k-space coordinates in a rosette MRSI 
%acquisition to account for any offsets due to gradient imperfections.  
%Each petal has N pts, and the central k-space point ideally coincides 
%with the first point in the adc, as well as all subsequent Nth points.  
%However, due to imperfections this is not the case.  
%What we need to do:  For each petal, we need to first estimate at which 
%point the maximum absolute signal intensity occurs.  This is going to be
%defined as our k=0 point for that petal.  Once we have done that, we are
%going to simply modify our k-file (leaving the data as-is) via a circular 
%shift, so that the k-coordinates now match the "true" locations in the a
%acquired dataset.  P.S. best to use the water unsuppressed data, as this
%provides the best SNR (also no need to combine averages).
%
% INPUTS:
% kFile    = The original k-file, which assumes k=0 at the first point of
%            the rosette trajectory 
% in       = The input MRSI dataset (preferrably water unsuppressed, prior
%            to any FT or density compensation).
%
% OUTPUTS:
% shiftedKTable = new k-table with the shifted coordinates
% ind           = vector of indices indicating the sample index of the 
%                 maximum k-space intensity (for each petal 


function [shiftedKTable, ind] = op_CSIrose_shiftKcoords(kFile, in)
    
%First, check if the input data has averages.  If so, combine them using
%op_CSIaverage:
if in.dims.averages
    in = op_CSIAverage(in);
end

%Now load the k-file:
[kTable,kArray] = readKFile(kFile);

%Find the number of k-space points per rosette cycle:
kPtsPerCycle = getKPtsPerCycle(kTable);

%Now loop through the coil elements and rosette petals to find the maximum
%signal intensity in a given cycle.  
for n=1:size(in.data,in.dims.coils)
    for m=1:size(in.data,in.dims.ky)
        [~,indMax(n,m)]=max(abs(in.data(1:kPtsPerCycle,n,m)));
    end
end

%Now take the average across all elements (we are using a single index for all TRs):
ind=round(mean(mean(indMax,1)));

%Now apply a circular shift to each TR in the kArray:
numTRs = size(kArray,1)/kPtsPerCycle;
kArrayShifted=zeros(size(kArray));
for n=1:numTRs
    tempX = kArray(((n-1)*kPtsPerCycle)+1:n*kPtsPerCycle,1);
    tempY = kArray(((n-1)*kPtsPerCycle)+1:n*kPtsPerCycle,2);
    kArrayShifted(((n-1)*kPtsPerCycle)+1:n*kPtsPerCycle,1) = circshift(tempX,ind-1);
    kArrayShifted(((n-1)*kPtsPerCycle)+1:n*kPtsPerCycle,2) = circshift(tempY,ind-1);
end

%Now put the shifted kArray back into a shifted kTable:
shiftedKTable = kTable;
shiftedKTable.Kx = kArrayShifted(:,1);
shiftedKTable.Ky = kArrayShifted(:,2);

%Now can we write this back into a k-file
writetable(shiftedKTable,[kFile(1:end-4) '_shifted.txt']);

