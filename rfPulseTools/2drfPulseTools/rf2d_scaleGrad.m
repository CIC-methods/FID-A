% rf2d_scaleGrad.m
% Jamie Near, 2025.
%
% USAGE:
% RF_out=rf2d_scaleGrad(2drf,scale);
% 
% DESCRIPTION:
% Scale the gradient amplitude of a gradient modulated 2DRF pulse by a fixed factor.  
% 
% INPUTS:
% 2drf      = Input 2DRF pulse definition structure.  Must be a gradient
%             modulated pulse.  
% scale     = Scaling factor by which you would like to multiply the 
%             existing gradient waveform.  
% 
% OUTPUTS:
% RF_out    = Output rf waveform following the scaling of gradient 
%             waveform.

function RF_out=rf2d_scaleGrad(rf2d,scale)

if ~isstruct(rf2d)
    error('ERROR:  the input 2DRF pulse must be in structure format.  Try using rf2d_create() to create it!  Aborting.  ');
end

newgrads = rf2d.gradients;
newgrads(:,1)=newgrads(:,1)*scale;
newgrads(:,2)=newgrads(:,2)*scale;

RF_out=rf2d;
RF_out.gradients=newgrads;
