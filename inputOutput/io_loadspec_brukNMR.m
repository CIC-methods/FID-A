%io_loadspec_brukNMR.m
%Chathura Kumaragamage, McGill University 2016.
%Jamie Near, McGill University 2016.
%
% USAGE:
% [out,ref]=io_loadspec_brukNMR(filename);
% 
% DESCRIPTION:
% Reads in Bruker MRS data (fid.raw, fid.ref).
%
% op_loadspec_bruk outputs the data in structure format, with fields corresponding to time
% scale, fids, frequency scale, spectra, and header fields containing
% information about the acquisition.  The resulting matlab structure can be
% operated on by the other functions in this MRS toolbox.
% 
% INPUTS:
% inDir   = Path to the scan number directory that contains the 'pdata' folder.
%
% OUTPUTS:
% out = Input dataset in FID-A structure format.
% ref = The Reference scan data (navigator echoes) in FID-A structure 
%       format, if applicable.



function [out,ref]=io_loadspec_brukNMR(inDir)


%get the original number of Repetitions (all repatitions from scanner is generated as
%a 1D vector. Need to split before further processing
averages=1; %Because Bruker does the averaging online.


%IMPORT USING FID FILE
fid_data=fread(fopen([inDir '/fid']),'int');
real_fid = fid_data(2:2:length(fid_data));
imag_fid = fid_data(1:2:length(fid_data));
fids=real_fid-1i*imag_fid;

%rawAverages=size(real_fid,1)./rawDataPoints;
rawAverages=1;
specs=fftshift(ifft(fids,[],1),1);
    
%calculate the size;
sz=size(specs);
out.flags.averaged=1;
%specify the dims
dims.t=1;
dims.coils=0;
dims.averages=0;
dims.subSpecs=0;
dims.extras=0;


%Now get some of the relevent spectral parameters

%First get the spectral width
acqu_fid=fopen([inDir '/acqu']);
line=fgets(acqu_fid);
index=findstr(line,'$SW_h=');
while isempty(index)
   line=fgets(acqu_fid);
   index=findstr(line,'$SW_h=');
end
equals_index=findstr(line,'=');
spectralwidth=line(equals_index+1:end);
spectralwidth=str2double(spectralwidth);
fclose(acqu_fid);

%Dwell time
dwelltime=1/spectralwidth;

%Now get the transmitter frequency
acqu_fid=fopen([inDir '/acqu']);
line=fgets(acqu_fid);
index=findstr(line,'$SFO1=');
while isempty(index)
   line=fgets(acqu_fid);
   index=findstr(line,'$SFO1=');
end
equals_index=findstr(line,'=');
txfrq=line(equals_index+1:end);
txfrq=str2double(txfrq);
txfrq=txfrq*1e6;
fclose(acqu_fid);

%Bo
Bo=txfrq/42577000;

%Spectral width in PPM
spectralwidthppm=spectralwidth/(txfrq/1e6);


%Now get the TE and TR
te=0;
tr=0;

%get the sequence
sequence=' ';


%Specify the number of subspecs.  For now, this will always be one.
subspecs=1;
rawSubspecs=1;


%calculate the ppm scale
ppm=[4.65+(spectralwidthppm/2):-spectralwidthppm/(length(specs)-1):4.65-(spectralwidthppm/2)];

%calculate the dwelltime:
dwelltime=1/spectralwidth;

%calculate the time scale
t=[0:dwelltime:(sz(1)-1)*dwelltime];



%FILLING IN DATA STRUCTURE FOR THE FID.RAW DATA
out.fids=fids;
out.specs=specs;
out.sz=sz;
out.ppm=ppm;  
out.t=t;    
out.spectralwidth=spectralwidth;
out.dwelltime=dwelltime;
out.txfrq=txfrq;
out.date=date;
out.dims=dims;
out.Bo=Bo;
out.averages=rawAverages;
out.rawAverages=rawAverages;
out.subspecs=subspecs;
out.rawSubspecs=rawSubspecs;
out.seq=sequence;
out.te=te;
out.tr=tr;
out.pointsToLeftshift=0;


%FILLING IN THE FLAGS FOR THE FID.RAW DATA
out.flags.writtentostruct=1;
out.flags.gotparams=1;
out.flags.leftshifted=0;
out.flags.filtered=0;
out.flags.zeropadded=0;
out.flags.freqcorrected=0;
out.flags.phasecorrected=0;

out.flags.addedrcvrs=1;
out.flags.subtracted=0;
out.flags.writtentotext=0;
out.flags.downsampled=0;
out.flags.avgNormalized=0;
if out.dims.subSpecs==0
    out.flags.isISIS=0;
else
    out.flags.isISIS=(out.sz(out.dims.subSpecs)==4);
end
