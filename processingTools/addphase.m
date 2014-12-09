%addphase.m
%Jamie Near, McGill University 2014.
%
%USAGE:
%PhasedSpecs=addphase(specs,AddedPhase);
%
%DESCRIPTION:  
%Add equal amount of complex phase to each point of a vector.
%
%INPUTS:
%specs          = Input vector.
%AddedPhase     = Amount of phase (degrees) to add.

function PhasedSpecs=addphase(specs,AddedPhase);

PhasedSpecs=specs.*(ones(size(specs))*exp(1i*AddedPhase*pi/180));