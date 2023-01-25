%rf_single2dualBand.m
%Peter Truong, Sunnybrook Research Institute 2021.
%
% USAGE:
% [rf,AMPINT]=rf_single2dualBand(rf_struct,tp,df);
% 
% DESCRIPTION:
% Creates a dualband RF pulse from user a single band RF pulse with a duration 
% of tp(ms). If the original rf pulse set at frequency = f_orig, the 
% second band will be at frequency = f_orig + df. 
% 
% INPUTS:
% rf_struct  = RF pulse definition structure.
% tp         = duration of the pulse from rf_struct in [ms].
% df         = frequency of 2nd band [Hz].
%
% OUTPUTS:
% rf         = Output rf waveform of the created dualband rf pulse, in FID-A rf 
%              pulse structure format.
% AMPINT     = Calculated amplitude integral (for use in Siemens .pta files).
function [rf,AMPINT]=rf_single2dualBand(rf_struct,tp,df)

% performing frequency shift on rf_struct
rf_struct2=rf_freqshift(rf_struct,tp,df);

% getting the waveforms in complex form
rf1_waveform=rf_struct.waveform(:,2).*exp(1i.*rf_struct.waveform(:,1)*pi/180);
rf2_waveform=rf_struct2.waveform(:,2).*exp(1i.*rf_struct2.waveform(:,1)*pi/180);

%Here we compute the Amplitude integral(AMPINT), which is used by magnetom
%to calculate the transmitter power that is required in order to achieve
%the desired flip angle
rf1_waveform_scaled=rf1_waveform./max(abs(rf1_waveform));
AI = sum(rf1_waveform_scaled);

combined_waveform=rf1_waveform+rf2_waveform;
combined_waveform_scaled=combined_waveform./max(abs(combined_waveform));
AMPINT=AI./max(abs(combined_waveform_scaled));

rf = rf_struct;
rf.waveform(:,1)=phase(combined_waveform_scaled).*180/pi;
rf.waveform(:,2)=abs(combined_waveform_scaled);

% updating the time-bandwith product & time-w1max product
rf.tbw=rf.tbw*2;
rf.tw1=rf.tw1*2;