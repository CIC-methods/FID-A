%io_writejmrui.m
%Jamie Near, McGill University 2014.
%Edits from
%   Edith Touchet-Valle, Texas A&M University 2024.
%   Jacob Degitz, Texas A&M University, 2024.
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
% outfile    = Desired filename of output text file. EV: Make sure to add '.txt' at the end of the filename
% Nuc        = Type of nucleus as a number. jMRUi offers 6 options:
%               0. Unknown
%               1. H
%               2. P
%               3. C
%               4. F
%               5. Na
% PatName    = Name of Patient
% a0         = zero order phase
% fpath      = *optional* data path
%
% OUTPUTS:
% RF         = Same as input.  Not used.  The primary output of this
%                function is a text file in jMRUI txt format. 

function RF=io_writejmrui(in,outfile,Nuc,PatName,zop,fpath);

Filename = in.fname;
t0=0;
Bo=in.Bo;
scanner='Verio';
addinfo='';

%type=input('Is this a difference spectrum (MEGA etc.)  y or n:  ','s');


% Determine data size
[num_data_pts,datasets] = size(in.specs);

% Force data path if not added explicitly
if nargin < 6 || isempty(fpath)
    fpath = dir(Filename).folder;
end




%write to txt file for jmrui
fid=fopen(fullfile([fpath,'\'],outfile),'w'); % Changed 5/24/24 by ETV
fprintf(fid,'jMRUI Data Textfile');
fprintf(fid,'\r\n');
fprintf(fid,'\r\nFilename: %s',Filename);
fprintf(fid,'\r\n');
fprintf(fid,'\r\nPointsInDataset: %i',num_data_pts);
fprintf(fid,'\r\nDatasetsInFile: %i',datasets);
fprintf(fid,'\r\nSamplingInterval: %4.6E',(1/in.spectralwidth)*1000);
fprintf(fid,'\r\nZeroOrderPhase: %1.0E',zop);
fprintf(fid,'\r\nBeginTime: %1.0E',t0);
fprintf(fid,'\r\nTransmitterFrequency: %4.6E',in.txfrq);
fprintf(fid,'\r\nMagneticField: %4.6E',Bo);
fprintf(fid,'\r\nTypeOfNucleus: %1.0E',Nuc);
fprintf(fid,'\r\nNameOfPatient: %s',PatName);
fprintf(fid,'\r\nDateOfExperiment: %s',strip(in.date_time(1:end-8)));
fprintf(fid,'\r\nSpectrometer: %s',scanner);
fprintf(fid,'\r\nAdditionalInfo: %s\r\n\r\n\r\n',addinfo);
fprintf(fid,'Signal\r\n');
fprintf(fid,'sig(real)\t sig(imag)\r\n');

fids=in.fids;
for n = 1:datasets
    fprintf(fid,'Signal %d out of %d in file\r\n',n,datasets);
    for row = 1:num_data_pts
        fprintf(fid,'%.8E\t %.8E\t %.8E\t %.8E\r\n', real(fids(row,n)), imag(fids(row,n)));
    end
end
fclose(fid);
RF = in;

