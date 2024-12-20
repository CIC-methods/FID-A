% addphase1.m
% Jamie Near, McGill University 2014.
% Edits from
%   Edith Touchet-Valle, Texas A&M University 2024.
%   Jacob Degitz, Texas A&M University 2024.
% 
% USAGE:
% PhasedSpecs=addphase1(specs,ppm,timeShift,B0,txfrq,ppm0);
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
% B0          = The main magnetic field strength (needed since ppm depends on
%               B0)
% txfrq       = The transmit frequency in MHz - added by ETV
% ppm0        = The frequency "origin" (ppm) of the 1st order phase shift. 
%               (this frequency will undergo 0 phase shift).
%
% OUTPUTS: 
% PhasedSpecs = Output vector (1st order phased version of the input). 

function PhasedSpecs=addphase1(specs,ppm,timeShift,B0,txfrq,ppm0);

if nargin<4
    gamma = (txfrq/1e6)/B0;
    if gamma > 42 % For proton
        ppm0 = 4.65;
    else % For all other nuclei
        ppm0 = 0;
    end
end


f=(ppm'-ppm0)*(txfrq); %Frequency scale in Hz
rep=size(specs);
rep(1)=1;
f=repmat(f,rep);

%we need to make a vector of added phase values based on the frequency
%scale of the spectrum.  the phase shift at any given point should be given
%by the time shift (in the time domain), multiplied by the frequency at
%that point:  phas[cycles]=f[Hz]*timeShift;
            % phas[radians]=f[Hz]*timeShift*2*pi;
phas=f*timeShift*2*pi;

PhasedSpecs=specs.*exp(-1i*phas);

%plot([1:b],PhasedSpecs(1,:));


