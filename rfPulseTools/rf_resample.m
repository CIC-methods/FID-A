% rf_resample.m
% Jamie Near, McGill University 2014.
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
newWaveform(:,1)=resample(RF_out.waveform(:,1),P,Q);
newWaveform(:,2)=resample(RF_out.waveform(:,2),P,Q);
newWaveform(:,3)=ones(length(newWaveform(:,1)),1);
if ~isPhsMod
    newWaveform(:,1)=180*(newWaveform(:,1)>100);
end
if RF_in.isGM 
    newWaveform(:,4)=resample(RF_out.waveform(:,4),P,Q);
end

RF_out.waveform=newWaveform;