% rf_timeReverse.m
% Jamie Near, Sunnybrook Research Institute 2022.
%
% USAGE:
% RF_out=rf_timeReverse(RF_in);
% 
% DESCRIPTION:
% Time-reverse the RF pulse.  
% 
% INPUTS:
% RF_in     = Input RF pulse definition structure.  
% 
% OUTPUTS:
% RF_out    = Time-reversed version of the input RF pulse.

function RF_out=rf_timeReverse(RF_in);

if ~isstruct(RF_in)
    error('ERROR:  the input RF pulse must be in structure format.  Try using rf_readwaveform to convert it!  Aborting.  ');
end

newWaveform=RF_in.waveform(end:-1:1,:);

RF_out=RF_in;
RF_out.waveform=newWaveform;
RF_out.rfCentre=1-RF_in.rfCentre;
