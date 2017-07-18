%io_writejmrui.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% RF=io_writejmrui(in,outfile);
% 
% DESCRIPTION:
% Takes MRS data in matlab structure format and writes it to a text file
% that can be read by jMRUI.
% 
% INPUTS:
% in         = input data in matlab structure format.
% outfile    = Desired filename of output text file.
%
% OUTPUTS:
% RF         = Same as input.  Not used.  The primary output of this
%                function is a text file in jMRUI txt format. 

function RF=io_writejmrui(in,outfile);

datsets=1;
zop=0;
t0=0;
Bo=in.Bo;
Nuc=0;
PatName='No Name';
scanner='TrioTim';
addinfo='jnear';

%type=input('Is this a difference spectrum (MEGA etc.)  y or n:  ','s');


%index=input('Enter Fid Index to use:  ');
RF=zeros(length(in.fids),4);
RF(:,1)=real(in.fids);
RF(:,2)=imag(in.fids);
RF(:,3)=real(in.specs);
RF(:,4)=imag(in.specs);




%write to txt file for jmrui
fid=fopen(outfile,'w+');
fprintf(fid,'jMRUI Data Textfile');
fprintf(fid,'\n\nFilename: %s' ,outfile);
fprintf(fid,'\n\nPointsInDataset: %i',length(RF(:,1)));
fprintf(fid,'\nDatasetsInFile: %i',datsets);
fprintf(fid,'\nSamplingInterval: %4.6E',in.dwelltime*1000);
fprintf(fid,'\nZeroOrderPhase: %1.0E',zop);
fprintf(fid,'\nBeginTime: %1.0E',t0);
fprintf(fid,'\nTransmitterFrequency: %4.6E',in.txfrq);
fprintf(fid,'\nMagneticField: %4.6E',Bo);
fprintf(fid,'\nTypeOfNucleus: %1.0E',Nuc);
fprintf(fid,'\nNameOfPatient: %s',PatName);
fprintf(fid,'\nDateOfExperiment: %i',in.date);
fprintf(fid,'\nSpectrometer: %s',scanner);
fprintf(fid,'\nAdditionalInfo: %s\n\n\n',addinfo);
fprintf(fid,'Signal and FFT\n');
fprintf(fid,'sig(real)\tsig(imag)\tfft(real)\tfft(imag)\n');
fprintf(fid,'Signal 1 out of %i in file\n',datsets);
fprintf(fid,'%1.8f\t%1.8f\t%1.8f\t%1.8f\n',RF');
fclose(fid);