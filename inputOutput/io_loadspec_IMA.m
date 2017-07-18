%io_loadspec_IMA.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% out=io_loadspec_IMA(filename,Bo,spectralwidth,te,tr);
% 
% DESCRIPTION:
% Loads a siemens .IMA file into matlab structure format.
% 
% INPUTS:
% filename       = Filename of Siemens .IMA file to load.
% Bo             = Field strength (Tesla).
% spectralwidth  = spectral width of the input spectrum (Hz).
% te             = Echo time (ms).  Optional.  Defulat is [].
% tr             = Repetition time (ms).  Optional.  Default is [].
%
% OUTPUTS:
% out        = Input dataset in FID-A structure format.

function out=io_loadspec_IMA(filename,Bo,spectralwidth,te,tr);

if nargin<5
    tr=[];
    if nargin<4
        te=[];
    end
end

%Load Dicom Info using Chris Rogers' "SiemensCsaParse.m" function:
info=SiemensCsaParse(filename);

%Read in Dicom file using Chris Rogers' "SiemensCsaReadFid.m" function:
[fids,info]=SiemensCsaReadFid(info,0);

sz=size(fids);

Naverages=1;
Ncoils=1;


if ndims(fids)==4  %Default config when 4 dims are acquired
    dims.t=1;
    dims.coils=2;
    dims.averages=3;
    dims.subSpecs=4;
elseif ndims(fids)<4  %To many permutations...ask user for dims.
    if Naverages == 1 && Ncoils == 1
        if ndims(fids)>1
            dims.t=1;
            dims.coils=0;
            dims.averages=0;
            dims.subSpecs=2;
        else
            dims.t=1;
            dims.coils=0;
            dims.averages=0;
            dims.subSpecs=0;
        end
    elseif Naverages>1 && Ncoils==1
        if ndims(fids)>2
            dims.t=1;
            dims.coils=0;
            dims.averages=2;
            dims.subSpecs=3;
        else
            dims.t=1;
            dims.coils=0;
            dims.averages=2;
            dims.subSpecs=0;
        end
    elseif Naverages==1 && Ncoils>1
        if ndims(fids)>2
            dims.t=1;
            dims.coils=2;
            dims.averages=0;
            dims.subSpecs=3;
        else
            dims.t=1;
            dims.coils=2;
            dims.averages=0;
            dims.subSpecs=0;
        end
    elseif Naverages>1 && Ncoils>1
        if ndims(fids)>3
            dims.t=1;
            dims.coils=2;
            dims.averages=3;
            dims.subSpecs=4;
        else
            dims.t=1;
            dims.coils=2;
            dims.averages=3;
            dims.subSpecs=0;
        end
    end
%     dims.t=1;
%     dims.coils=input('Enter the coils Dimension (0 for none):  ');
%     dims.averages=input('Enter the averages Dimension (0 for none):  ');
%     dims.subSpecs=input('Enter the subSpecs Dimension (0 for none);  ');
end

specs=fftshift(ifft(fids,[],dims.t),dims.t);


    

%Now get relevant scan parameters:*****************************

%Calculate Dwell Time
dwelltime=1/spectralwidth;

%Calculate TxFrq
txfrq=42577000*Bo;

%Get Date
date=01012000;

%Find the number of averages.  'averages' will specify the current number
%of averages in the dataset as it is processed, which may be subject to
%change.  'rawAverages' will specify the original number of acquired 
%averages in the dataset, which is unchangeable.
if dims.subSpecs ~=0
    if dims.averages~=0
        averages=sz(dims.averages)*sz(dims.subSpecs);
        rawAverages=averages;
    else
        averages=sz(dims.subSpecs);
        rawAverages=1;
    end
else
    if dims.averages~=0
        averages=sz(dims.averages);
        rawAverages=averages;
    else
        averages=1;
        rawAverages=1;
    end
end

%Find the number of subspecs.  'subspecs' will specify the current number
%of subspectra in the dataset as it is processed, which may be subject to
%change.  'rawSubspecs' will specify the original number of acquired 
%subspectra in the dataset, which is unchangeable.
if dims.subSpecs ~=0
    subspecs=sz(dims.subSpecs);
    rawSubspecs=subspecs;
else
    subspecs=1;
    rawSubspecs=subspecs;
end

%****************************************************************


%Calculate t and ppm arrays using the calculated parameters:
f=[(-spectralwidth/2)+(spectralwidth/(2*sz(1))):spectralwidth/(sz(1)):(spectralwidth/2)-(spectralwidth/(2*sz(1)))];
ppm=f/(Bo*42.577);
ppm=ppm+4.65;

t=[0:dwelltime:(sz(1)-1)*dwelltime];


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
out.seq='';
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
if out.dims.subSpecs==0
    out.flags.isISIS=0;
else
    out.flags.isISIS=(out.sz(out.dims.subSpecs)==4);
end



%DONE
