% getKPtsPerCycle.m
% Brenden Kadota, Sunnybrook Research Institute 2021.
% Jamie Near, Sunnybrook Research Institute 2024.
%
% USAGE: 
% kPtsPerCycle = getKPtsPerCycle(kTable)
%
% DESCRIPTION:
%This function accepts as input a "kTable" (the output of the function 
%'readkFile.m' for a rapid MRSI seuqence and returns the number of k-space 
%points sampled in a single cycle of the trajectory.  A cycle is defined as 
%the k-space pattern that is traced repeatedly in a single TR of a MRSI 
%sequence.  For example, for concentric rings or rosette MRSI, one cycle 
%would be a single concentric ring or rosette petal, respectively, and the 
%output numKPointsPerCycle would be the number of ADC samples collected in 
%one pass around the ring/petal.  In the case of EPSI, the output 
%numKPointsPerCycle would be the number of adc samples collected between 
%repeated passes of the acquired line.  
%
%Note:  This function was previously called "calculateNumSpatialPoints.m",
%but the name was changed to reflect the fact that these are not spatial
%points being collected, but rather k-space samples.  
%
%INPUTS:
% kTable                 = Matlab table specifying the k-space coordinates 
%                          of a rapid MRSI acquisition.  This kTable can be
%                          obtained by reading a k-space coordinates text
%                          file using the function 'readKFile.m'.  
%
%OUTPUTS
% kPtsPerCycle:   = The calculated number of k-space points per cycle.


function kPtsPerCycle = getKPtsPerCycle(kTable)
        
    %Find the maximum TR index in the trajectory:
    num_TR = max(kTable.TR, [], 'all');
    
    %Divide the number of lines in the table by the number of TRs to get
    %the number of points per TR/cycle:
    kPtsPerCycle = height(kTable)/num_TR;

end