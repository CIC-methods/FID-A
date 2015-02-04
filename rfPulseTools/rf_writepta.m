%rf_writepta.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% RF=rf_writepta(rf,outfile);
% 
% DESCRIPTION:
% Write a matlab RF pulse structure (containing 3 x N waveform array field 
% with rf.waveform(1,:)= phase, rf.waveform(2,:)=amplitude, and rf.waveform(3,:)=timestep), 
% to a siemens format .pta file.
% 
% INPUTS:
% rf         = matlab RF pulse.
% outfile    = name of the output .pta file to be written.

function RF=rf_writepta(rf,outfile);

B1INT=sum(rf.waveform(:,2));

l=length(rf.waveform(:,2));
RF.waveform=zeros(3,l);
RF.waveform(1,:)=rf.waveform(:,2)/max(rf.waveform(:,2));
RF.waveform(2,:)=mod(rf.waveform(:,1)*pi/180,2*pi);
RF.waveform(3,:)=[0:l-1];


%I copied this from the code in svs_se.cpp... I think it's just integrating
%under the absolute part of the pulse.  
POWERINT=sum(RF.waveform(1,:).*(cos(RF.waveform(2,:))+sin(RF.waveform(2,:)))); 

%AMPLITUDE INTEGRAL (AMPINT)...  This is the parameter that is used to
%calculate the transmitter voltage at runtime, and therefore it determines
%the flip angle that will be realied.  It must be calculated properly to
%ensure that the flip angle will be correct.  In addition, the calculation
%of AMPINT depends on the type of pulse that you are dealing with.  For a
%straigt plain wave (That is, the phase of each data point is either a
%constant ph, or ph+180 or ph-180)  the AMPINT is simply the sum of the
%waveform points when the maximum value is scaled to 1.  For example, for a
%three point pulse with point values of -1, 3 and -1, the AMPINT would be
%1.  This will be the method used for automatic calculation of AMPINT.  If
%the pulse is not a plane pulse, do not use this method to calculate
%AMPINT.

AIcalc=input('Would you like to perform an automatic Amplitude Inegral Calculation?(y or n)  ', 's');

if AIcalc=='y' || AIcalc=='Y'
    R=sum(RF.waveform(1,:).*cos(RF.waveform(2,:)))
    I=sum(RF.waveform(1,:).*sin(RF.waveform(2,:)))
    AMPINT=sqrt((R^2)+(I^2))
end


%If the pulse is not a plane wave (Adiabatic pulses, off resonance pulses, 
%or dual Band Pulses), you should calculate the Amplitude integral 
%yourself and input it manually.  For methods of calculating AMPINT in this 
%case, see the document entitled 
%"RF_amplitude_integral_in_pulse_sequences.pdf', which can be found in
%/Users/jnear/Documents/RF pulses/ 

if AIcalc=='n' || AIcalc=='N'
    AMPINT=input('Now you must input the Amplitude Integral Manually.  Lets have it?  ');
end


%Reference gradient is only used in slice selective pulses.  This refers to
%the gradient required to excite a 10mm slice if the pulse duration is
%5.12ms  Given that the Time-bandwidth product is included in the rf
%structure, the refgrad value can be calculated here, but the option is
%provided to enter it manually, in case a different refgrad value is
%desired.  
refgrad=1;
refg=input('Would you like to manually input a reference gradient?  \nThis is the gradient amplitude in [mT/m] required to excite a 10mm slice with a \npulse duration of 5.12 ms.  If you say ''n'' (recommended for conventional pulses), \nthen this program will calculate refgrad for you. (y or n)  ','s');
if refg=='y' || refg=='Y'
    refgrad=input('Enter the reference gradient (mT/m):  ');
elseif refg=='n' || refg=='N'
    %refgrad [mT/m] = (1000 [mT]/[1T]) * tbw [unitless] ) / ( (0.00512 [s]) * (42577000 [Hz]/[T]) * (0.01 [m]) );
    refgrad = 1000 * RF.tbw/(0.00512 * 42577000 * 0.01);
end


%write to pta file for siemens
fid=fopen(outfile,'w+');
fprintf(fid,'REFGRAD: \t%5.6f' ,refgrad);
fprintf(fid,'\nAMPINT: \t%5.6f',AMPINT);
fprintf(fid,'\nPOWERINT: \t%5.6f\n',POWERINT);
n=0;
fprintf(fid,'%1.8f\t%1.8f\t;(%2.0f)\n',RF.waveform);
fclose(fid);