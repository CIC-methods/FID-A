% rf_addGrad.m
% Jamie Near, McGill University 2020.
%
% USAGE:
% RF_out=rf_addGrad(RF_in,grad);
% 
% DESCRIPTION:
% Add a gradient waveform to the input RF pulse. 
% 
% INPUTS:
% RF_in     = Input RF pulse definition structure
% grad      = If grad is a scalar, then rf_addGrad will add a constant
%             gradient with amplitude equal to the value of grad.  If grad 
%             is a vector, then rf_addGrad will add the specified gradient 
%             vector to the RF waveform.  In either case, the gradient 
%             should be specified in units of [G/cm].
% 
% OUTPUTS:
% RF_out    = Output rf waveform following the addition of gradient 
%             waveform.

function RF_out=rf_addGrad(RF_in,grad);

if ~isstruct(RF_in)
    error('ERROR:  the input RF pulse must be in structure format.  Try using rf_readwaveform to convert it!  Aborting.  ');
end

newWaveform=RF_in.waveform;

%Check if the input RF pulse already has a gradient. 
if size(newWaveform,2)>3
    keepGoing='n';
    disp('WARNING:  Input waveform already has a gradient waveform!');
    keepGoing=input('Do you wish to overwrite the existing gradient waveform (y or [n])','s');
    
    if strcmp(keepGoing,'y') || strcmp(keepGoing,'Y');
        disp('OK.  Overwriting existing gradient waveform.');
    elseif strcmp(keepGoing,'n') || strcmp(keepGoing,'N');
        error('ABORTING!!!');
    else
        error('Response not recognized.  ABORTING!!');
    end
end

%Now add the gradient. 
if isscalar(grad)
    %Adding a scalar gradient
    newWaveform(:,4)=grad*ones(size(newWaveform,1),1);
else
    %Adding a gradient waveform
    if size(grad) ~= size(newWaveform,1)
        if size(grad') ~= size(newWaveform,1);
            error('ERROR:  Gradient waveform does not match the length of the input RF pulse waveform!!  ABORTING!!');
        end
        newWaveform(:,4)=grad';
    end
    newWaveform(:,4)=grad';
end

RF_out=RF_in;
RF_out.waveform=newWaveform;