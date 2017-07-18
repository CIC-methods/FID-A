%io_writeRF.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% RF=io_writeRF(rf,outfile);
% 
% DESCRIPTION:
% Write a matlab RF pulse structure (containing 3 x N waveform array field 
% with rf.waveform(1,:)= phase, rf.waveform(2,:)=amplitude, and rf.waveform(3,:)=timestep), 
% to a varian/agilent format .RF file.
% 
% INPUTS:
% rf         = matlab RF pulse.
% outfile    = name of the output .RF file to be written.
%
% OUTPUTS:
% RF         = Same as input.  Not used.  The primary output of this
%                function is a text file in Varian/Agilent .RF format. 

function RF=io_writeRF(rf,outfile);

AM=rf.waveform(:,2);
ph=rf.waveform(:,1);
N=length(AM);
B1INT=sum(AM);

l=length(rf.waveform(:,2));
RF.waveform=zeros(3,l);
RF.waveform(1,:)=rf.waveform(:,2)/max(rf.waveform(:,2));
RF.waveform(2,:)=rf.waveform(:,1)*pi/180;
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

% AIcalc=input('Would you like to perform an automatic Amplitude Inegral Calculation?(y or n)  ', 's');
% 
% if AIcalc=='y' || AIcalc=='Y'
%     R=sum(RF(1,:).*cos(RF(2,:)))
%     I=sum(RF(1,:).*sin(RF(2,:)))
%     AMPINT=sqrt((R^2)+(I^2))
% end


%If the pulse is not a plane wave (Adiabatic pulses, off resonance pulses, 
%or dual Band Pulses), you should calculate the Amplitude integral 
%yourself and input it manually.  For methods of calculating AMPINT in this 
%case, see the document entitled 
%"RF_amplitude_integral_in_pulse_sequences.pdf', which can be found in
%/Users/jnear/Documents/RF pulses/ 

% if AIcalc=='n' || AIcalc=='N'
%     AMPINT=input('Now you must input the Amplitude Integral Manually.  Lets have it?  ');
% end

% Generate RF pulse data
    [x,npts]=size(RF.waveform);
    integral = sum(AM)/(max(AM)*npts);
    pulse = abs(AM);
    pulse=[ph; 1023*pulse./max(pulse)];
    

    %outfile=input('pick a name for the .RF file  ','s');
    
    if ~isempty(outfile)
        if length(outfile) < 3
             % Append .RF ending, outfile too short to have it
            outfile = [outfile '.RF'];
        elseif outfile(end-2:end) == '.rf'
            % Change .rf ending to .RF ending
            outfile(end-2:end) = '.RF';
        elseif outfile(end-2:end) ~= '.RF'
            % Append .RF if it does not exist
            outfile = [outfile '.RF'];
        end
	
        % Check if file exists and confirm overwrite
        if exist(outfile, 'file')
            r = input(['Delete existing file ' outfile ' (Y/n)?'],'s');
            if ~isempty(r)
                if lower(r(1)) == 'n'
                    return
                end
            end
        end
   
        fid=fopen(outfile,'w+');
        if fid < 0
            error('SLR: Unable to create %s.\n', outfile);
        end
    end
    
    %R=(R1+R2)/2;
    if ~isempty(outfile)
        % Write header to RF file
        fprintf(fid,'# RF Pulse Created by Jamie Near\n');
        %fprintf(fid,'# R= %5.2f\n',R);
        fprintf(fid,'# number of points =     %5.2f\n',N);
        fprintf(fid,'# TYPE        selective\n');
        fprintf(fid,'# MODULATION  amplitude\n');
        %fprintf(fid,'# EXCITEWIDTH %5.2f\n',R);
        %fprintf(fid,'# INVERTWIDTH %5.2f\n',R);
        fprintf(fid,'# INTEGRAL    -1\n');
        fprintf(fid,'# ******************************************************\n');
	
        fprintf(fid,'%8.3f  %9.3f  1.0\n',pulse);
        fprintf(fid,'#\n');
        fclose(fid);
    end