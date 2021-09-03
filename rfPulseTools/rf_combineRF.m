%rf_combineRF.m
%Peter Truong, Sunnybrook Research Institute 2021.
%
% USAGE:
% [rf,AMPINT]=rf_combineRF(rf_struct1,rf_struct2);
% 
% DESCRIPTION:
% Combines two similar rf waveforms (e.g. same pulse, but one frequency
% shifted), in FID-A rf pulse structure format.
% 
% INPUTS:
% rf_struct1 = RF pulse definition structure.
% rf_struct2 = RF pulse definition structure, similar to rf_struct1

%
% OUTPUTS:
% rf         = Output rf waveform of the combined rf pulse, in FID-A rf 
%              pulse structure format.
% AMPINT     = Calculated amplitude integral (for use in Siemens .pta files).
function [rf,AMPINT]=rf_combineRF(rf_struct1,rf_struct2)

%Convert waveforms to complex form
rf1_waveform=rf_struct1.waveform(:,2).*exp(1i.*rf_struct1.waveform(:,1)*pi/180);
rf2_waveform=rf_struct2.waveform(:,2).*exp(1i.*rf_struct2.waveform(:,1)*pi/180);

%Here we compute the Amplitude integral(AMPINT), which is used by magnetom
%to calculate the transmitter power that is required in order to achieve
%the desired flip angle
rf1_waveform_scaled=rf1_waveform./max(abs(rf1_waveform));
AI=sum(rf1_waveform_scaled);

combined_waveform=rf1_waveform + rf2_waveform;
combined_waveform_scaled=combined_waveform./max(abs(combined_waveform));
AMPINT=AI./max(abs(combined_waveform_scaled));

rf=rf_struct1;
rf.waveform(:,1)=phase(combined_waveform_scaled).*180/pi;
rf.waveform(:,2)=abs(combined_waveform_scaled);

% updating time-bandwidth product & time-w1max product
rf.tbw=rf_struct1.tbw+rf_struct2.tbw;
rf.tw1=rf_struct1.tw1+rf_struct2.tw1;
