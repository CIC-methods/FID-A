% addphase1.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% PhasedSpecs=addphase1(specs,ppm,timeShift,ppm0,B0);
% 
% DESCRIPTION:
% Add first order phase to a spectrum (added phase is linearly dependent on
% frequency).  This function operates on a vector (fid or spectrum), not on 
% a FID-A data structure.  For a phase shifting function that operates on a 
% FID-A data structure, see 'op_addphase.m'.
% 
% INPUTS:
% specs       = input vector
% ppm         = freqeuncy scale (ppm) corresponding to the specs vector
% timeShift   = This defines the amount of 1st order phase shift by
%               specifying the equivalent horizontal shift (in seconds) in the time
%               domain.
% ppm0        = The frequency "origin" (ppm) of the 1st order phase shift. 
%               (this frequency will undergo 0 phase shift).
% B0          = The main magnetic field strength (needed since ppm depends on
%               B0)
%
% OUTPUTS: 
% PhasedSpecs = Output vector (1st order phased version of the input). 

function PhasedSpecs=addphase1(specs,ppm,timeShift,ppm0,B0);

if nargin<4
    ppm0=4.65;
end


f=(ppm'-ppm0)*42.577*B0; %Frequency scale in Hz;  Assumes proton.
rep=size(specs);
rep(1)=1;
f=repmat(f,rep);

%we need to make a vector of added phase values based on the frequency
%scale of the spectrum.  the phase shift at any given point should be given
%by the time shift (in the time domain), multiplied by the frequency at
%that point:  phas[cycles]=f[Hz]*timeShift;
            % phas[radians]=f[Hz]*timeShift*2*pi;
phas=f*timeShift*2*pi;

PhasedSpecs=specs.*exp(-i*phas);

%plot([1:b],PhasedSpecs(1,:));


