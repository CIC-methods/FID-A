%io_writeRFtxt.m
%Jamie Near, McGill University 2020.
%
% USAGE:
% RF=io_writeRFtxt(rf,outfile);
% 
% DESCRIPTION:
% Write a matlab RF pulse structure (containing N X 3 (or 4) waveform array 
%field with rf.waveform(:,1)= phase [deg], rf.waveform(:,2)=amplitude 
%[a.u.], and rf.waveform(:,3)=timestep), and optionally rf.waveform(:,4)= 
%gradient [G/cm], to an RF pulse text file with 4 columns in the following 
%order:  Amplitude, phase, time-step, and gradient.
%  
% INPUTS:
% rf         = matlab RF pulse.
% outfile    = name of the output .RF file to be written.
%
% OUTPUTS:
% RF         = RF pulse array with columns in order of RF amplitude, phase, 
%              time-step and gradient strength.  This output is only for 
%              verification purposes.  The primary output of this function 
%              is a text file with columns formated in the same order.

function RF=io_writeRFtxt(rf,outfile);

RFtemp=rf.waveform;
RF=RFtemp;
RF(:,1)=RFtemp(:,2);
RF(:,2)=RFtemp(:,1);

N=length(RF(:,1));

% Generate RF pulse data
        
fid=fopen(outfile,'w+');
if fid < 0
    error('ERROR: Unable to create %s.\n', outfile);
end
    
if ~isempty(outfile)
    for n=1:N
        fprintf(fid,'%5.5f  %5.5f  %5.5f  %5.5f\n',RF(n,:));
    end
%    fprintf(fid,'#\n');
    fclose(fid);
end