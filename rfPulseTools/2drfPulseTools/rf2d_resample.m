% rf2d_resample.m
% Jamie Near, 2025.
% tl-edit: edited by Thomas Lange, University of Freiburg
%
% USAGE:
% RF_out=rf2d_resample(RF_in,N);
% 
% DESCRIPTION:
% Resample the input RF pulse into a new waveform with N discrete points.  
% 
% INPUTS:
% rf2d     = Input RF pulse definition structure
% N         = Number of points in new RF waveform
% 
% OUTPUTS:
% RF_out    = Output rf waveform following resampling.

function RF_out=rf2d_resample(rf2d,N)

if ~isstruct(rf2d)
    error('ERROR:  the input RF pulse must be in structure format. Aborting.  ');
end

P=N;
Q=length(rf2d.waveform(:,2));

%Find out if the pulse is phase modulated.  Whether or not it is phase
%modulated will determine how we go about resampling the phase function
a=(round(rf2d.waveform(:,2))==0)|(round(rf2d.waveform(:,2))==6*pi);

if sum(a)<length(rf2d.waveform(:,2))
    isPhsMod=true;
else
    isPhsMod=false;
end

RF_out=rf2d;

% tl-edit++ 
% convert waveform to complex array for resampling
waveform_compl=RF_out.waveform(:,1).*exp(1i*RF_out.waveform(:,2));

waveform_compl_new=resample(waveform_compl,P,Q);

newWaveform(:,2)=round(10000*angle(waveform_compl_new))/10000 + pi;
newWaveform(:,1)=(round(10000*abs(waveform_compl_new))/10000) / max(round(10000*abs(waveform_compl_new))/10000);


if ~isPhsMod
    newWaveform(:,2)=round(newWaveform(:,2));
    %newWaveform(:,2)=newWaveform(:,2)+(newWaveform(:,2) == -1*pi)*2*pi;
end

newGradients(:,1)=resample(RF_out.gradients(:,1),P,Q);
newGradients(:,2)=resample(RF_out.gradients(:,2),P,Q);

RF_out.waveform=newWaveform;
RF_out.gradients = newGradients;