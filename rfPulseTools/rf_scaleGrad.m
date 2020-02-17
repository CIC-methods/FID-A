% rf_scaleGrad.m
% Jamie Near, McGill University 2020.
%
% USAGE:
% RF_out=rf_scaleGrad(RF_in,scale);
% 
% DESCRIPTION:
% Scale the gradient amplitude of a gradient modulated pulse by a fixed factor.  
% 
% INPUTS:
% RF_in     = Input RF pulse definition structure.  Must be a gradient
%             modulated pulse.  
% scale     = Scaling factor by which you would like to multiply the 
%             existing gradient waveform.  
% 
% OUTPUTS:
% RF_out    = Output rf waveform following the scaling of gradient 
%             waveform.

function RF_out=rf_scaleGrad(RF_in,scale);

if ~isstruct(RF_in)
    error('ERROR:  the input RF pulse must be in structure format.  Try using rf_readwaveform to convert it!  Aborting.  ');
end

if ~RF_in.isGM
    error('ERROR:  the input RF pulse must be a gradient modulated pulse.  ABORTING!');
end

newWaveform=RF_in.waveform;
newWaveform(:,4)=newWaveform(:,4)*scale;

RF_out=RF_in;
RF_out.waveform=newWaveform;
RF_out.tthk=RF_in.tthk/scale;
