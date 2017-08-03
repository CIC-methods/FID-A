%io_loadspec_bruk.m
%Chathura Kumaragamage, McGill University 2016.
%Jamie Near, McGill University 2016.
%
% USAGE:
% [out,ref]=io_loadspec_bruk(filename);
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



function [out,ref]=io_loadspec_bruk(inDir)

%%%%%%%%chathu mod starts

%get the original number of Repetitions (all repatitions from scanner is generated as
%a 1D vector. Need to split before further processing
averages=1; %Because Bruker does the averaging online.
ADC_OFFSET=70;      %(offset between ADC on and ADC acquire, empirically determined)
method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'$PVM_DigNp=');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'$PVM_DigNp=');
end
equals_index=findstr(line,'=');
rawDataPoints=line(equals_index+1:end);
rawDataPoints=str2double(rawDataPoints);
fclose(method_fid);

%NOW LOAD IN THE RAW DATA.  FIRST TRY USING THE FID.RAW FILE.  IF THAT DOES
%NOT WORK, THEN USE THE REGULAR FID FID.
try
    fid_data=fread(fopen([inDir '/fid.raw']),'int');
    real_fid = fid_data(2:2:length(fid_data));
    imag_fid = fid_data(1:2:length(fid_data));
    fids_raw=real_fid-1i*imag_fid;
    
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
    fids_raw=reshape(fids_raw,[],rawAverages);
    fids_trunc=fids_raw(ADC_OFFSET:end,:);
    %         fids_trunc=fids_trunc';
    fids=padarray(fids_trunc, [ADC_OFFSET-1,0],'post');
    specs=fftshift(ifft(fids,[],1),1);
    sz=size(specs);
    out.flags.averaged=0;
    %specify the dims
    dims.t=1;
    dims.coils=0;
    dims.averages=2;
    dims.subSpecs=0;
    dims.extras=0;
catch
    %NOW TRY USING FID FILE
    display('WARNING: /fid.raw not found. Using /fid ....')
    %specs=fread(fopen([inDir '/pdata/1/2dseq']),'int');
    fid_data=fread(fopen([inDir '/fid']),'int');
    real_fid = fid_data(2:2:length(fid_data));
    imag_fid = fid_data(1:2:length(fid_data));
    fids_raw=real_fid-1i*imag_fid;
    
    rawAverages=size(real_fid,1)./rawDataPoints;
    
    if mod(size(real_fid,1),rawDataPoints) ~=0
        display 'number of repetitions cannot be accurately found';
    end
    fids_raw=reshape(fids_raw,rawDataPoints,[]);
    
    fids_trunc=fids_raw(ADC_OFFSET:end,:);
    fids=padarray(fids_trunc, [ADC_OFFSET-1,0],'post');
    %convert back to time domain
    %if the length of Fids is odd, then you have to do a circshift of one to
    %make sure that you don't introduce a small frequency shift into the fids
    %vector.
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
end

%NOW TRY LOADING IN THE REFERENCE SCAN DATA (IF IT EXISTS)
isRef=false;
if exist([inDir '/fid.ref'])
    isRef=true;
    ref_data=fread(fopen([inDir '/fid.ref']),'int');
    real_ref = ref_data(2:2:length(ref_data));
    imag_ref = ref_data(1:2:length(ref_data));
    fids_ref=real_ref-1i*imag_ref;
    
    fids_ref=reshape(fids_ref,[],rawAverages);
    fids_ref_trunc=fids_ref(ADC_OFFSET:end,:);
    %         fids_trunc=fids_trunc';
    reffids=padarray(fids_ref_trunc, [ADC_OFFSET-1,0],'post');
    refspecs=fftshift(ifft(reffids,[],1),1);
    sz=size(refspecs);
    out.flags.averaged=0;
    %specify the dims
    refdims.t=1;
    refdims.coils=0;
    refdims.averages=2;
    refdims.subSpecs=0;
    refdims.extras=0;
end


%%%%%Chathu mod ends


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



%Now get the TE
method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'$PVM_EchoTime=');
%%%%%CHATHU MOD for FIDS with Pulse acquire
loop_count=0;
TE_present=true;
while isempty(index) && TE_present
    line=fgets(method_fid);
    index=findstr(line,'$PVM_EchoTime=');
    loop_count=loop_count+1;
    if loop_count>10000
        TE_present=false;
    end
end

if TE_present
    equals_index=findstr(line,'=');
    te=line(equals_index+1:end);
    te=str2double(te);
else
    te=0;
end

%%%%%CHATHU MOD for FIDS with Pulse acquire END
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


if isRef
    %FILLING IN DATA STRUCTURE FOR THE FID.REF DATA
    ref.fids=reffids;
    ref.specs=refspecs;
    ref.sz=sz;
    ref.ppm=ppm;
    ref.t=t;
    ref.spectralwidth=spectralwidth;
    ref.dwelltime=dwelltime;
    ref.txfrq=txfrq;
    ref.date=date;
    ref.dims=dims;
    ref.Bo=Bo;
    ref.averages=rawAverages;
    ref.rawAverages=rawAverages;
    ref.subspecs=subspecs;
    ref.rawSubspecs=rawSubspecs;
    ref.seq=sequence;
    ref.te=te;
    ref.tr=tr;
    ref.pointsToLeftshift=0;
    
    
    %FILLING IN THE FLAGS FOR THE FID.REF DATA
    ref.flags.writtentostruct=1;
    ref.flags.gotparams=1;
    ref.flags.leftshifted=0;
    ref.flags.filtered=0;
    ref.flags.zeropadded=0;
    ref.flags.freqcorrected=0;
    ref.flags.phasecorrected=0;
    
    ref.flags.addedrcvrs=1;
    ref.flags.subtracted=0;
    ref.flags.writtentotext=0;
    ref.flags.downsampled=0;
    ref.flags.avgNormalized=0;
    if ref.dims.subSpecs==0
        ref.flags.isISIS=0;
    else
        ref.flags.isISIS=(ref.sz(ref.dims.subSpecs)==4);
    end
end
