% rf_resample.m
% Jamie Near, McGill University 2014.
% tl-edit: edited by Thomas Lange, University of Freiburg
%
% USAGE:
% RF_out=rf_resample(RF_in,N);
% 
% DESCRIPTION:
% Resample the input RF pulse into a new waveform with N discrete points.  
% 
% INPUTS:
% RF_in     = Input RF pulse definition structure
% N         = Number of points in new RF waveform
% 
% OUTPUTS:
% RF_out    = Output rf waveform following resampling.

function RF_out=rf_resample(RF_in,N);

if ~isstruct(RF_in)
    error('ERROR:  the input RF pulse must be in structure format.  Try using rf_readwaveform to convert it!  Aborting.  ');
end

P=N;
Q=length(RF_in.waveform(:,1));

%Find out if the pulse is phase modulated.  Whether or not it is phase
%modulated will determine how we go about resampling the phase function
a=(round(RF_in.waveform(:,1))==180)|(round(RF_in.waveform(:,1))==0);

if sum(a)<length(RF_in.waveform(:,1))
    isPhsMod=true;
else
    isPhsMod=false;
end

RF_out=RF_in;

% tl-edit++ 
% convert waveform to complex array for resampling
waveform_compl=RF_out.waveform(:,2).*exp(i*RF_out.waveform(:,1)*pi/180);

waveform_compl_new=resample(waveform_compl,P,Q);

newWaveform(:,1)=round(10000*angle(waveform_compl_new)*180/pi)/10000;
newWaveform(:,2)=round(10000*abs(waveform_compl_new))/10000;
newWaveform(:,3)=ones(length(newWaveform(:,1)),1);

if ~isPhsMod
    newWaveform(:,1)=round(newWaveform(:,1));
    newWaveform(:,1)=newWaveform(:,1)+(newWaveform(:,1) == -180)*360;
end
% tl-edit--

if RF_in.isGM 
    newWaveform(:,4)=resample(RF_out.waveform(:,4),P,Q);
end

RF_out.waveform=newWaveform;