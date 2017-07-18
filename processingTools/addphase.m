% addphase.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% PhasedSpecs=addphase(specs,AddedPhase);
% 
% DESCRIPTION:  
% Add equal amount of complex phase to each point of a vector.  This
% function operates on a vector (fid or spectrum), not on a FID-A data
% structure.  For a phase shifting function that operates on a FID-A data
% structure, see 'op_addphase.m'.
% 
% INPUTS:
% specs          = Input vector.
% AddedPhase     = Amount of phase (degrees) to add.
%
% OUTPUTS: 
% PhasedSpecs    = Output vector (0th order phased version of the input). 


function PhasedSpecs=addphase(specs,AddedPhase);

PhasedSpecs=specs.*(ones(size(specs))*exp(1i*AddedPhase*pi/180));