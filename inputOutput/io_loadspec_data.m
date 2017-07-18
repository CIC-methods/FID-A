%io_loadspec_data.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% [out,out_w]=io_loadspec_data(filename,sw,Larmor,subspecs,te,tr);
% 
% DESCRIPTION:
% Reads in philips MRS data (.data and .list files) using code adapted from 
% PhilipsRead_data.m, provided as part of the Gannet software package by 
% Richard Edden (gabamrs.blogspot.com).
% 
% op_loadspec_data outputs the data in structure format, with fields 
% corresponding to time scale, fids, frequency scale, spectra, and header 
% fields containing information about the acquisition.  The resulting 
% matlab structure can be operated on by the other functions in this MRS 
% toolbox.  NOTE:  Since the Gannet code is geared towards edited GABA MRS 
% data, this code may not be general enough to handle all types of MRS data.  
% Suggestions are most welcome.
% 
% INPUTS:
% filename   = filename of Philips .data file to be loaded.
% sw         = spectral width (Hz) 
% Larmor     = Larmor frequency (Hz/ppm, ie.  127 for 3T)
% subspecs   = number of subspectra in the data (from spectral editing, ISIS, etc.)
% te         = echo time (ms).  Optional, default is [].
% tr         = repetition time (ms).  Optional, default is [].
%
% OUTPUTS:
% out        = Input water suppressed dataset in FID-A structure format.
% out_w      = Input water reference dataset in FID-A structure format. 


function [out,out_w]=io_loadspec_data(filename,sw,Larmor,subspecs,te,tr);

if nargin<6
    if nargin<5
        te=[];
    end
    tr=[];
end


%read in the data using the philipsDataLoad.m (adapted from PhilipsRead_data.m)
[FullData,WaterData]=philipsDataLoad(filename);

%As far as I can tell, the data that comes out of the philipsDataLoad
%function is normally a N x Navgs matrix.  The coils have already been 
%combined. The Navgs dimension contains all the subspectra, so we will 
%split them now.  Note, that in the data-list format that I have seen, the
%edit-OFF subspectra appear in elements [1 2 5 6 9 10 13 14...] and the
%edit-ON subspectra appear in the elements [3 4 7 8 11 12 15 16...].  Other sequences
%may result in a different subspecs order, but for now we will separate the 
%subspectra in this way.
%If the data has multiple subspectra 
if subspecs>1
    %First make an vector that holds the indices of the ON subspectra:
    totalAvgs=size(FullData,2);
    OFFindices=[1:2:totalAvgs]-mod([0:(totalAvgs/2)-1],2);
    ONindices=[2:2:totalAvgs]+mod([1:totalAvgs/2],2);
    %Now split the subspectra out of the "averages" dimension:
    data(:,:,1)=FullData(:,OFFindices);
    data(:,:,2)=FullData(:,ONindices);
else
    data=FullData;
end

fids=squeeze(data);
fids_w=squeeze(WaterData)';

sz=size(fids);
sz_w=size(fids_w);

%Find the magnetic field strength:
Bo=Larmor/42.577;

%Find the number of averages:
Naverages=size(fids,2)*size(fids,3);
Naverages_w=size(fids_w,2)*size(fids_w,3);

%In Philips data/list format, coil channels have already been combined:
Ncoils=1;
Ncoils_w=1;

%Now create a record of the dimensions of the data array.  
dims.t=1;
dims.coils=0;
dims.averages=2;
if subspecs>1
    dims.subSpecs=3;
else
    dims.subSpecs=0;
end

dims_w.t=1;
dims_w.coils=0;
dims_w.averages=2;
dims_w.subSpecs=0;


specs=fftshift(ifft(fids,[],dims.t),dims.t);
specs_w=fftshift(ifft(fids_w,[],dims_w.t),dims_w.t);


%Now get relevant scan parameters:*****************************

%Get Spectral width and Dwell Time
spectralwidth=sw;
dwelltime=1/spectralwidth;
    
%Get TxFrq
txfrq=Larmor*1e6;

%Leave date blank
date='';

%Find the number of averages.  'averages' will specify the current number
%of averages in the dataset as it is processed, which may be subject to
%change.  'rawAverages' will specify the original number of acquired 
%averages in the dataset, which is unchangeable.
%FOR WATER SUPPRESSED DATA:
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

%FOR WATER UNSUPPRESSED DATA:
if dims_w.subSpecs ~=0
    if dims_w.averages~=0
        averages_w=sz(dims_w.averages)*sz(dims_w.subSpecs);
        rawAverages_w=averages_w;
    else
        averages_w=sz(dims_w.subSpecs);
        rawAverages_w=1;
    end
else
    if dims_w.averages~=0
        averages_w=sz(dims_w.averages);
        rawAverages_w=averages_w;
    else
        averages_w=1;
        rawAverages_w=1;
    end
end


%Find the number of subspecs.  'subspecs' will specify the current number
%of subspectra in the dataset as it is processed, which may be subject to
%change.  'rawSubspecs' will specify the original number of acquired 
%subspectra in the dataset, which is unchangeable.
%FOR WATER SUPPRESSED DATA:
if dims.subSpecs ~=0
    subspecs=sz(dims.subSpecs);
    rawSubspecs=subspecs;
else
    subspecs=1;
    rawSubspecs=subspecs;
end

%FOR WATER UNSUPPRESSED DATA:
if dims_w.subSpecs ~=0
    subspecs_w=sz(dims.subSpecs);
    rawSubspecs_w=subspecs_w;
else
    subspecs_w=1;
    rawSubspecs_w=subspecs_w;
end

%****************************************************************


%Calculate t and ppm arrays using the calculated parameters:
f=[(-spectralwidth/2)+(spectralwidth/(2*sz(1))):spectralwidth/(sz(1)):(spectralwidth/2)-(spectralwidth/(2*sz(1)))];
ppm=f/(Bo*42.577);
ppm=ppm+4.65;

t=[0:dwelltime:(sz(1)-1)*dwelltime];


%FOR WATER SUPPRESSED DATA
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
out.flags.averaged=0;
out.flags.addedrcvrs=0;
out.flags.subtracted=0;
out.flags.writtentotext=0;
out.flags.downsampled=0;
if out.dims.subSpecs==0
    out.flags.isISIS=0;
else
    out.flags.isISIS=(out.sz(out.dims.subSpecs)==4);
end


%FOR WATER UNSUPPRESSED DATA
%FILLING IN DATA STRUCTURE
out_w.fids=fids_w;
out_w.specs=specs_w;
out_w.sz=sz_w;
out_w.ppm=ppm;  
out_w.t=t;    
out_w.spectralwidth=spectralwidth;
out_w.dwelltime=dwelltime;
out_w.txfrq=txfrq;
out_w.date=date;
out_w.dims=dims_w;
out_w.Bo=Bo;
out_w.averages=averages_w;
out_w.rawAverages=rawAverages_w;
out_w.subspecs=subspecs_w;
out_w.rawSubspecs=rawSubspecs_w;
out_w.seq='';
out_w.te=te;
out_w.tr=tr;
out_w.pointsToLeftshift=0;


%FILLING IN THE FLAGS
out_w.flags.writtentostruct=1;
out_w.flags.gotparams=1;
out_w.flags.leftshifted=0;
out_w.flags.filtered=0;
out_w.flags.zeropadded=0;
out_w.flags.freqcorrected=0;
out_w.flags.phasecorrected=0;
out_w.flags.averaged=0;
out_w.flags.addedrcvrs=0;
out_w.flags.subtracted=0;
out_w.flags.writtentotext=0;
out_w.flags.downsampled=0;
if out_w.dims.subSpecs==0
    out_w.flags.isISIS=0;
else
    out_w.flags.isISIS=(out.sz(out.dims.subSpecs)==4);
end



%DONE
