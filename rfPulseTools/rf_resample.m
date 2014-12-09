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
% RF        = RF pulse definition structure
% N         = Number of points in new RF waveform


function RF_out=rf_resample(RF_in,N);

if ~isstruct(RF_in)
    error('ERROR:  the input RF pulse must be in structure format.  Try using rf_readwaveform to convert it!  Aborting.  ');
end

P=1;
Q=ceil(length(RF_in.waveform(:,1))/N);


RF_out=RF_in;
newWaveform(:,1)=resample(RF_out.waveform(:,1),P,Q);
newWaveform(:,2)=resample(RF_out.waveform(:,2),P,Q);
newWaveform(:,3)=ones(length(newWaveform(:,1)),1);
newWaveform(:,1)=180*(newWaveform(:,1)>100);

RF_out.waveform=newWaveform;