%io_loadspec_irBruk.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% out=io_loadspec_irBruk(inDir);
% 
% DESCRIPTION:
% Reads in Bruker MRS data (1i and 1r files).  Generally with this data,
% the averages and coil channels have already been combined.
%
% op_loadspec_irBruk outputs the data in structure format, with fields corresponding to time
% scale, fids, frequency scale, spectra, and header fields containing
% information about the acquisition.  The resulting matlab structure can be
% operated on by the other functions in this MRS toolbox.
% 
% INPUTS:
% inDir   = Path to the scan directory that contains the 'pdata' folder.
%
% OUTPUTS:
% out        = Input dataset in FID-A structure format.

function out=io_loadspec_irBruk(inDir)

real_fileID=fopen([inDir '/pdata/1/1r']);
imag_fileID=fopen([inDir '/pdata/1/1i']);

real=fread(real_fileID,'int');
imag=fread(imag_fileID,'int');

specs=real-1i*imag;

%convert back to time domain
%if the length of Fids is odd, then you have to do a circshift of one to
%make sure that you don't introduce a small frequency shift into the fids
%vector.
if mod(length(specs),2)==0
    %disp('Length of vector is even.  Doing normal conversion');
    fids=fft(fftshift(specs,1),[],1);
else
    %disp('Length of vector is odd.  Doing circshift by 1');
    fids=fft(circshift(fftshift(specs,1),1),[],1);
end


%calculate the size;
sz=size(fids);

%Now get some of the relevent spectral parameters

%First get the spectral width
method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'$PVM_DigSw=');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'$PVM_DigSw=');
end
equals_index=findstr(line,'=');
spectralwidth=line(equals_index+1:end);
spectralwidth=str2double(spectralwidth);
fclose(method_fid);

%Now get the transmitter frequency
acqp_fid=fopen([inDir '/acqp']);
line=fgets(acqp_fid);
index=findstr(line,'$BF1=');
while isempty(index)
    line=fgets(acqp_fid);
    index=findstr(line,'$BF1=');
end
equals_index=findstr(line,'=');
txfrq=line(equals_index+1:end);
txfrq=str2double(txfrq);
txfrq=txfrq*1e6;
fclose(acqp_fid);

%B0
Bo=txfrq/42577000;

%Spectral width in PPM
spectralwidthppm=spectralwidth/(txfrq/1e6);

%Now get the original number of averages
averages=1; %Because Bruker does the averaging online.
method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'$PVM_NAverages=');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'$PVM_NAverages=');
end
equals_index=findstr(line,'=');
rawAverages=line(equals_index+1:end);
rawAverages=str2double(rawAverages);
fclose(method_fid);

%Now get the TE
method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'$PVM_EchoTime=');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'$PVM_EchoTime=');
end
equals_index=findstr(line,'=');
te=line(equals_index+1:end);
te=str2double(te);
fclose(method_fid);

%Now get the TR
method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'$PVM_RepetitionTime=');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'$PVM_RepetitionTime=');
end
equals_index=findstr(line,'=');
tr=line(equals_index+1:end);
tr=str2double(tr);
fclose(method_fid);

%Now get the sequence
method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'$Method=');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'$Method=');
end
equals_index=findstr(line,'=');
sequence=line(equals_index+1:end);
sequence=strtrim(sequence);
fclose(method_fid);

%Specify the number of subspecs.  For now, this will always be zero.
subspecs=0;
rawSubspecs=0;


%calculate the ppm scale
ppm=[4.65+(spectralwidthppm/2):-spectralwidthppm/(length(specs)-1):4.65-(spectralwidthppm/2)];

%calculate the dwelltime:
dwelltime=1/spectralwidth;

%calculate the time scale
t=[0:dwelltime:(sz(1)-1)*dwelltime];

%specify the dims
dims.t=1;
dims.coils=0;
dims.averages=0;
dims.subSpecs=0;
dims.extras=0;

%FILLING IN DATA STRUCTURE
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
out.averages=averages;
out.rawAverages=rawAverages;
out.subspecs=subspecs;
out.rawSubspecs=rawSubspecs;
out.seq=sequence;
out.te=te;
out.tr=tr;
out.pointsToLeftshift=0;


%FILLING IN THE FLAGS
out.flags.writtentostruct=1;
out.flags.gotparams=1;
out.flags.leftshifted=0;
out.flags.filtered=0;
out.flags.zeropadded=0;
out.flags.freqcorrected=0;
out.flags.phasecorrected=0;
out.flags.averaged=1;
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

