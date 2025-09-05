% rf_2dtimeReverse.m
% Jamie Near, Sunnybrook Research Institute 2022.
%
% USAGE:
% RF_out=rf_2dtimeReverse(rf2d);
% 
% DESCRIPTION:
% Time-reverse the 2DRF pulse.  
% 
% INPUTS:
% rf2d    = Input RF pulse definition structure.  
% 
% OUTPUTS:
% RF_out    = Time-reversed version of the input RF pulse.
%
% NOTES: Pulses put into this function CANNOT have a ramp! Remove ramp in rf2d_create() through
% setting 'withRamp = False'.

function RF_out = rf2d_timeReverse(rf2d)

    if ~isstruct(rf2d)
        error('ERROR: input RF pulse must be a struct. Use rf2d_create() to generate the waveform.');
    end

    amp = rf2d.waveform(:,1);

    % Count leading and trailing zeros only
    leadingZeros  = find(amp,1,'first') - 1;
    counter = leadingZeros;

    if counter > 2
        error(['ERROR!! Cannot reverse a ramped pulse. ' ...
               'Set withRamp = False in rf2d_create() and try remaking the pulse again.']);
    else
        RF_out = rf2d;
        RF_out.waveform  = rf2d.waveform(end:-1:1,:);
        RF_out.gradients = rf2d.gradients(end:-1:1,:);
    end
end
