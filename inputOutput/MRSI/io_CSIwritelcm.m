%io_writelcm.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% RF=io_CSIwritelcm(in,outfile,te);
% 
% DESCRIPTION:
% Takes MRS data in matlab structure format and writes it to a text file
% that can be read by LCModel.
% 
% INPUTS:
% in         = input data in matlab structure format.
% outfile    = Desired filename of output text file.
% te         = Echo time of acquisition (in ms).
%
% OUTPUTS:
% RF         = Same as input.  Not used.  The primary output of this
%                function is a text file in LCModel raw format. 

function RF=io_CSIwritelcm(in,outfile)


if ~in.flags.addedrcvrs
    error('ERROR:  reciever channels must be combined first');
end


% datsets=1;
% zop=0;
% t0=0;
  Bo=in.Bo;
  hzppm=getGamma('overTwoPi', true)*Bo/1e6;
  dwelltime=in.spectralDwellTime;
% Nuc=0;
% PatName='No Name';
% scanner='TrioTim';
% addinfo='jnear';

specs = ifft(fftshift(getData(in), in.dims.t), [], in.dims.t);
vec_signal = specs(:);

RF = [imag(vec_signal) real(vec_signal)];


%write to txt file for jmrui
fid=fopen(outfile,'w+');
fprintf(fid,' $SEQPAR');
%fprintf(fid,'\n\nFilename: %s' ,outfile);
%fprintf(fid,'\n\nPointsInDataset: %i',length(data_struct.fids));
%fprintf(fid,'\nDatasetsInFile: %i',datsets);
%fprintf(fid,'\nSamplingInterval: %1.0E',data_struct.dwelltime*1000);
%fprintf(fid,'\nZeroOrderPhase: %1.0E',zop);
%fprintf(fid,'\nBeginTime: %1.0E',t0);
%fprintf(fid,'\nTransmitterFrequency: %2.4E',data_struct.txfrq);
%fprintf(fid,'\nMagneticField: %2.1E',Bo);
%fprintf(fid,'\nTypeOfNucleus: %1.0E',Nuc);
%fprintf(fid,'\nNameOfPatient: %s',PatName);
%fprintf(fid,'\nDateOfExperiment: %i',data_struct.date);
%fprintf(fid,'\nSpectrometer: %s',scanner);
%fprintf(fid,'\nAdditionalInfo: %s\n\n\n',addinfo);
%fprintf(fid,'Signal and FFT\n');
%fprintf(fid,'sig(real)\tsig(imag)\tfft(real)\tfft(imag)\n');
%fprintf(fid,'Signal 1 out of %i in file\n',datsets);
fprintf(fid,'\n echot= %2.2f',in.te);
fprintf(fid,'\n seq= ''PRESS''');
fprintf(fid,'\n hzpppm= %5.6f',hzppm);
fprintf(fid,'\n NumberOfPoints= %i',in.sz(1));
fprintf(fid,'\n dwellTime= %5.6f' ,dwelltime);
fprintf(fid,'\n $END');
fprintf(fid,'\n $NMID');
fprintf(fid,'\n id=''ANONYMOUS '', fmtdat=''(2E15.6)''');
fprintf(fid,'\n volume=8.0');
fprintf(fid,'\n tramp=1.0');
fprintf(fid,'\n $END\n');
fprintf(fid,'  % 7.6e  % 7.6e\n',RF');
fclose(fid);
