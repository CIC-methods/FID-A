%io_loadspec_bruk.m
%Jamie Near, Sunnybrook Research Institute, 2023
%Wendy Oakden, Sunnybrook Research Institute, 2020
%Chathura Kumaragamage, McGill University 2016.
%
% USAGE:
% [out,ref]=io_loadspec_bruk(inDir, rawData, );
% 
% DESCRIPTION:
% Reads in Bruker MRS data (fid.raw, fid.ref).
%
% op_loadspec_bruk reads a Bruker data file and outputs a data structure in 
% FID-A format (with fields corresponding to time scale, fids, frequency 
% scale, spectra, and header fields containing information about the 
% acquisition.  The resulting matlab structure can be operated on by the 
% other FID-A functions.
%
% This function gives you the option to read in either the raw data (with
% individual averages stored separately), or the combined data. 
%
% For Bruker versions PV5 and earlier, the water suppressed and water
% unsuppressed data were acquired separately.  In such cases, 
% io_loadspec_bruk would need to be run separately on the water suppressed 
% and water unsuppressed scan folders, and the 2nd output 'ref' refers to 
% any separately acquired navigator echoes.  For Bruker versions PV6 and later, 
% it became common for the water suppressed data to be acquired
% concurrently with the water unsuppressed data (within a single scan ID).
% For such datasets, this function reads both the water suppressed data
% (1st output, "out"), and the water unsuppressed data (2nd output
% "ref").
% 
% INPUTS:
% inDir     = String variable specifying the path to the scan number 
%             directory.  Alternatively, inDir can be an integer specifying 
%             number of the scan directory to analyze (assuming it is in the
%             current directory).  
% rawData   = 'y' or 'n' (default = 'y')
%             'y' - load the individually stored fids 
%             'n' - load the Bruker combined data
% leftshift = Number of points to left-shift the FID.  (default = 68, the
%             number of points normally added/acquired before the echo on 
%             Bruker systems.) 
%
% OUTPUTS:
% out = Input dataset in FID-A structure format.
% ref = The Reference scan data (Navigator echoes in PV5, water reference 
%       data in later versions) in FID-A structure format, if applicable.

function [out,ref]=io_loadspec_bruk(inDir, rawData, leftshift)

if nargin<3
    leftshift=68; %default number of points before the echo in most Bruker Data
    if nargin<2
        rawData='y';
        if nargin<1
            error('ERROR: no input spectrum specified.  Aborting!!');
        end
    end
end

%Allow the user to pass the input directory as either a string or an
%integer.
if isnumeric(inDir)
    if ~rem(inDir,1)==0
        error('ERROR: only integer values of ''inDir'' are allowed.  ABORTING');
    else
        inDir=num2str(inDir);
    end
end

%Find out what version of ParaVision we are working with:
acqp_fid=fopen([inDir '/acqp']);
line=fgets(acqp_fid);
index=strfind(line,'<PV');
while isempty(index)
    line=fgets(acqp_fid);
    index=strfind(line,'<PV');
end
version=line(2:end-2);

%get the number of raw data points
method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
if contains(version,'PV-360')
    index=strfind(line,'$PVM_SpecMatrix=');
    while isempty(index)
        line=fgets(method_fid);
        index=strfind(line,'$PVM_SpecMatrix=');
    end
    line=fgets(method_fid);
    rawDataPoints=str2double(line);
else
    index=strfind(line,'$PVM_DigNp=');
    while isempty(index)
        line=fgets(method_fid);
        index=strfind(line,'$PVM_DigNp=');
    end
    equals_index=strfind(line,'=');
    rawDataPoints=line(equals_index+1:end);
    rawDataPoints=str2double(rawDataPoints);
end
fclose(method_fid);

%Find the number of averages in the main dataset (number of averages in
%the 'ref' dataset will be determined later):
method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=strfind(line,'$PVM_NAverages=');
while isempty(index)
    line=fgets(method_fid);
    index=strfind(line,'$PVM_NAverages=');
end
equals_index=strfind(line,'=');
rawAverages=line(equals_index+1:end);
rawAverages=str2double(rawAverages);
fclose(method_fid);

%Find if multiple RF channels are used in the main dataset
method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=strfind(line,'$PVM_EncUseMultiRec=');
while isempty(index)
    line=fgets(method_fid);
    index=strfind(line,'$PVM_EncUseMultiRec=');
end
equals_index=strfind(line,'=');
multiRec=line(equals_index+1:end-1);
if strcmp(multiRec,'Yes')
    multiRcvrs=true;
else
    multiRcvrs=false;
end
%Now determine the number of RF channels used:
index=strfind(line,'$PVM_EncNReceivers=');
while isempty(index)
    line=fgets(method_fid);
    index=strfind(line,'$PVM_EncNReceivers=');
end
equals_index=strfind(line,'=');
Nrcvrs=line(equals_index+1:end);
Nrcvrs=str2num(Nrcvrs);
fclose(method_fid);

%Now load in the data.  Either raw or averaged data, depending on what was
%requested via the 'rawData' parameter. 
if strcmp(rawData,'y') || strcmp(rawData,'Y')
    if contains(version,'PV 6') || contains(version,'PV 7') || contains(version,'PV-360')
        data = fopen([inDir '/rawdata.job0']); % WO - changed for PV6.0
        fid_data = fread(data,'int32');
    elseif contains(version,'PV 5')
        data = fopen([inDir '/fid.raw']);
        fid_data = fread(data,'int');
    end
    real_fid = fid_data(2:2:length(fid_data));
    imag_fid = fid_data(1:2:length(fid_data));
    fids_raw=real_fid-1i*imag_fid;
    fclose(data); %WO - added fclose

    averages=rawAverages;  %since these data are uncombined;
    out.flags.averaged=0; %make the flags structure

    %Reshape to put the averages along a 2nd dimension
    fids_raw=reshape(fids_raw,[],rawAverages);

    %If there are multiple receivers *I think* that these always get stored
    %separately by default in the fid.raw file.  Therefore, at this stage, 
    %if this is a PV360 dataset with multiple receivers, we need to reshape 
    %the dataset again:
    if ~contains(version,'PV 5') && multiRcvrs
        fids_raw=reshape(fids_raw,rawDataPoints,Nrcvrs,rawAverages);
        %Permute so that time is along 1st dimension, averages is along 2nd 
        %dimension, and coils is along 3rd dimension:
        fids_raw=permute(fids_raw,[1,3,2]); 
    end
elseif strcmp(rawdata,'n') || strcmp(rawData,'N')
    %REQUEST PROCESSED DATA ONLY:  Use the FID file instead of fid.raw.
    data = fopen([inDir '/fid']);
    fid_data=fread(data,'int');
    real_fid = fid_data(2:2:length(fid_data));
    imag_fid = fid_data(1:2:length(fid_data));
    fids_raw=real_fid-1i*imag_fid;
    fclose(data); %WO

    averages=1; %since we requested the combined data.
    out.flags.averaged=1; %make the flags structure


    %***JN***
    % REMOVED THIS SECTION.  DON'T THINK IT IS NECESSARY.
    %     if mod(size(real_fid,1),rawDataPoints) ~=0
    %         display 'number of repetitions cannot be accurately found';
    %     end
    %     fids_raw=reshape(fids_raw,rawDataPoints,[]);
    %***JN***
else
    error('ERROR:  rawData variable not recognized.  Options are ''y'' or ''n''.');
end

%Perform left-shifting to remove points before the echo
fids_trunc=fids_raw(leftshift+1:end,:,:);

%replace the left-shifted points with zeros at the end
fids=padarray(fids_trunc, [leftshift,0],'post');

%Do the fourier transform
%specs=fftshift(ifft(fids,[],1),1);
specs=FIDAfft(fids,1,'t');
sz=size(specs); %size of the array

%specify the dims:
%Time dimension:
dims.t=1;%the time dimension is always the 1st dimension

%Coils dimension:
%As far as I am aware, the RF coils are only stored separately in PV360.
if ~contains(version,'PV 5') && multiRcvrs
    %Coils dimension should normally be after the averages dimension,
    %unless there are no averages, in which case the coild dimension will
    %be after the time dimension.  
    if rawAverages==1
        dims.coils=2;
    elseif rawAverages>1
        dims.coils=3;
    end
else
    dims.coils=0;
end

if strcmp(rawData,'y') || strcmp(rawData,'Y')
    dims.averages=2;
elseif strcmp(rawData,'n') || strcmp(rawData,'N')
    dims.averages=0;
end
%I have not encountered any Bruker datasets so far where there are
%subSpectra or extras dimensions.  
dims.subSpecs=0;
dims.extras=0;


%NOW TRY LOADING IN THE REFERENCE SCAN DATA (IF IT EXISTS)
isRef=false;
if contains(version,'PV 5')
    isRef=exist([inDir '/fid.ref']);
elseif contains(version,'PV 6') || contains(version,'PV 7') || contains(version,'PV-360.2')
    isRef=exist([inDir '/fid.refscan']);
elseif contains(version,'PV-360.3')
    isRef=exist([inDir '/pdata/1/fid_refscan.64']);
end

if isRef
    if contains(version,'PV 5')
        data = fopen([inDir '/fid.ref']);
        ref_data=fread(data,'int16');
    elseif contains(version,'PV 6') || contains(version,'PV 7') || contains(version,'PV-360.2')
        data = fopen([inDir '/fid.refscan']);
        ref_data=fread(data,'int32');
    elseif contains(version,'PV-360.3')
        data = fopen([inDir '/pdata/1/fid_refscan.64']);
        ref_data=fread(data,'float64');
    end
    real_ref = ref_data(2:2:length(ref_data));
    imag_ref = ref_data(1:2:length(ref_data));
    fids_ref=real_ref-1i*imag_ref;
    fclose(data); %WO - added fclose

    %Find the number of averages in the ref dataset:
    if contains(version,'PV 5')
        rawAverages_ref=rawAverages;
        averages_ref=rawAverages_ref;
    else
        method_fid=fopen([inDir '/method']);
        line=fgets(method_fid);
        index=strfind(line,'$PVM_RefScanNA=');
        while isempty(index)
            line=fgets(method_fid);
            index=strfind(line,'$PVM_RefScanNA=');
        end
        equals_index=strfind(line,'=');
        rawAverages_ref=line(equals_index+1:end);
        rawAverages_ref=str2double(rawAverages_ref);
        fclose(method_fid);
        averages_ref=rawAverages_ref;
    end

    fids_ref=reshape(fids_ref,[],rawAverages_ref);
    fids_ref_trunc=fids_ref(leftshift+1:end,:);
    reffids=padarray(fids_ref_trunc, [leftshift,0],'post');
    %refspecs=fftshift(ifft(reffids,[],1),1);
    refspecs=FIDAfft(reffids,1,'t');
    sz_ref=size(refspecs);
    ref.flags.averaged=0;
    %specify the dims
    refdims.t=1;
    refdims.coils=0;
    if rawAverages_ref>1
        refdims.averages=2;
    else
        refdims.averages=0;
    end    
    refdims.subSpecs=0;
    refdims.extras=0;
else
    %Ref scans not found.  Print warning:
    disp('WARNING REFERENCE SCANS NOT FOUND.  RETURNING EMPTY REF STRUCTURE.');
end

%Now get some of the relevent spectral parameters
%First get the spectral width
method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
if contains(version,'PV 5') || contains(version,'PV 6') || contains(version,'PV 7')
    index=strfind(line,'$PVM_DigSw=');
    while isempty(index)
        line=fgets(method_fid);
        index=strfind(line,'$PVM_DigSw=');
    end
    equals_index=strfind(line,'=');
    spectralwidth=line(equals_index+1:end);
    spectralwidth=str2double(spectralwidth);
elseif contains(version,'PV-360')
    index=strfind(line,'$PVM_SpecSWH=');
    while isempty(index)
        line=fgets(method_fid);
        index=strfind(line,'$PVM_SpecSWH=');
    end
    line=fgets(method_fid);
    spectralwidth=str2double(line);

    %finding nucleus tag
    index_nuc=findstr(line,'$PVM_Nucleus1Enum=');
    if index_nuc
        line_nuc=line;
    end
end

% getting the nucleus and gamma
equals_index=strfind(line_nuc,'=');
nucleus=line_nuc(equals_index+1:end-1);
gamma=getgamma(nucleus);

fclose(method_fid);


%Now get the transmitter frequency
acqp_fid=fopen([inDir '/acqp']);
line=fgets(acqp_fid);
index=strfind(line,'$BF1=');
while isempty(index)
    line=fgets(acqp_fid);
    index=strfind(line,'$BF1=');
end
equals_index=strfind(line,'=');
txfrq=line(equals_index+1:end);
txfrq=str2double(txfrq);
txfrq=txfrq*1e6;
fclose(acqp_fid);

%B0
Bo=txfrq/(gamma*1e6);

% %Spectral width in PPM
% spectralwidthppm=spectralwidth/(txfrq/1e6);

%Now get the TE
method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=strfind(line,'$PVM_EchoTime=');
while isempty(index) && ~feof(method_fid)
    line=fgets(method_fid);
    index=strfind(line,'$PVM_EchoTime=');
end
if ~feof(method_fid)
    equals_index=strfind(line,'=');
    te=line(equals_index+1:end);
    te=str2double(te);
else
    te=0;
end
fclose(method_fid);

%Now get the TR
method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=strfind(line,'$PVM_RepetitionTime=');
while isempty(index)
    line=fgets(method_fid);
    index=strfind(line,'$PVM_RepetitionTime=');
end
equals_index=strfind(line,'=');
tr=line(equals_index+1:end);
tr=str2double(tr);
fclose(method_fid);

%Now get the sequence
method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=strfind(line,'$Method=');
while isempty(index)
    line=fgets(method_fid);
    index=strfind(line,'$Method=');
end
equals_index=strfind(line,'=');
sequence=line(equals_index+1:end);
sequence=strtrim(sequence);
fclose(method_fid);

%Specify the number of subspecs.  For now, this will always be one.
subspecs=1;
rawSubspecs=1;

%calculate the ppm scale
% ppm=[4.65+(spectralwidthppm/2):-spectralwidthppm/(length(specs)-1):4.65-(spectralwidthppm/2)];
ppm=calcppm(spectralwidth,sz(1),Bo,gamma);
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
out.averages=averages;
out.rawAverages=rawAverages;
out.subspecs=subspecs;
out.rawSubspecs=rawSubspecs;
out.seq=sequence;
out.te=te;
out.tr=tr;
out.pointsToLeftshift=68-leftshift;
out.version=version;

%PT - 2023
out.nucleus=nucleus;
out.gamma=gamma;
%Unsure how to set header and filename, update later - PT,2023
out.hdr='';
out.filename=inDir;

%FILLING IN THE FLAGS FOR THE FID.RAW DATA
out.flags.writtentostruct=1;
out.flags.gotparams=1;
out.flags.leftshifted=0;
out.flags.filtered=0;
out.flags.zeropadded=0;
out.flags.freqcorrected=0;
out.flags.phasecorrected=0;

if multiRcvrs && dims.coils
    out.flags.addedrcvrs=0;
else
    out.flags.addedrcvrs=1;
end
out.flags.subtracted=0;
out.flags.writtentotext=0;
out.flags.downsampled=0;
out.flags.avgNormalized=0;
if out.dims.subSpecs==0
    out.flags.isFourSteps=0;
else
    out.flags.isFourSteps=(out.sz(out.dims.subSpecs)==4);
end


if isRef
    %FILLING IN DATA STRUCTURE FOR THE FID.REF DATA
    ref.fids=reffids;
    ref.specs=refspecs;
    ref.sz=sz_ref;
    ref.ppm=ppm;
    ref.t=t;
    ref.spectralwidth=spectralwidth;
    ref.dwelltime=dwelltime;
    ref.txfrq=txfrq;
    ref.date=date;
    ref.dims=refdims;
    ref.Bo=Bo;
    ref.averages=averages_ref;
    ref.rawAverages=rawAverages_ref;
    ref.subspecs=subspecs;
    ref.rawSubspecs=rawSubspecs;
    ref.seq=sequence;
    ref.te=te;
    ref.tr=tr;
    ref.pointsToLeftshift=68-leftshift;
    ref.version=version;
    
    %PT - 2023
    ref.nucleus=nucleus;
    ref.gamma=gamma;
    %Unsure how to set header and filename, update later - PT,2023
    ref.hdr='';
    ref.filename=inDir;
    
    
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
        ref.flags.isFourSteps=0;
    else
        ref.flags.isFourSteps=(ref.sz(ref.dims.subSpecs)==4);
    end
else
    %REF NOT FOUND.  RETURN EMPTY STRUCTURE.
    ref = [];
end