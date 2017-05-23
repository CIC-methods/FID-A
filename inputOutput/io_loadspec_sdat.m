%io_loadspec_sdat.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% out=io_loadspec_sdat(filename,subspecs);
% 
% DESCRIPTION:
% Reads in Philpis MRS data (.spar and .sdat files) using code adapted from 
% PhilipsRead.m, provided as part of the Gannet software package by Richard 
% Edden (gabamrs.blogspot.com).
% 
% op_loadspec_sdat outputs the data in structure format, with fields 
% corresponding to time scale, fids, frequency scale, spectra, and header 
% fields containing information about the acquisition.  The resulting 
% matlab structure can be operated on by the other functions in this MRS 
% toolbox.  NOTE:  Since the Gannet code is geared towards edited GABA MRS 
% data, this code may not be general enough to handle all types of MRS 
% data.  Suggestions are most welcome.  ALSO:  This code is not currently
% smart enough to parse out all of the relevent information from the header
% file, such as the number of subspectra.  So for now, these details must 
% be passed to the function as input arguments.  Help implementing these 
% improvements are most welcome!!
% 
% INPUTS:
% filename   = filename of Philips sdat file to be loaded.
% subspecs   = number of subspectra in the data (from spectral editing, ISIS, etc.)

function out=io_loadspec_sdat(filename,subspecs);

%read in the data using the GELoad.m (adapted from GERead.m)
philipsOut=philipsLoad(filename);

%As far as I can tell, the data that comes out of the philipsLoad
%function is normally a N x Navgs x Ncoils matrix.  The Navgs dimension
%contains all the subspectra, so we will split them now:
%If the data has multiple subspectra 
if subspecs==2
    %Split the subspectra out of the "averages" dimension:
    data(:,:,1)=philipsOut(:,[1:2:end-1],:);
    data(:,:,2)=philipsOut(:,[2:2:end],:);
elseif subspecs==4
    data(:,:,1)=philipsOut(:,[1:4:end-3],:);
    data(:,:,2)=philipsOut(:,[2:4:end-2],:);
    data(:,:,3)=philipsOut(:,[3:4:end-1],:);
    data(:,:,4)=philipsOut(:,[4:4:end],:);
else
    data=philipsOut;
end

fids=squeeze(data);

sz=size(fids);

%Read in the spar file to get important parameters:
sparheader=textread([filename(1:end-4) 'spar'],'%s');

%Find the centre frequency
txfrq_index=find(ismember(sparheader, 'synthesizer_frequency')==1);
txfrq = str2num(sparheader{txfrq_index+2});

%Calculate what the magnetic field strength was:
Bo=txfrq/42577000;

%Find the number of averages:
Naverages=size(fids,2)*size(fids,3);

%Find the spectral width
spectralwidth_index=find(ismember(sparheader, 'sample_frequency')==1);
spectralwidth = str2num(sparheader{spectralwidth_index+2});

%Calculate the dwell time
dwelltime=1/spectralwidth;

%Find the echo time
te_index=find(ismember(sparheader, 'echo_time')==1);
te=str2num(sparheader{te_index+2});

%Find the repetition time
tr_index=find(ismember(sparheader, 'repetition_time')==1);
tr=str2num(sparheader{tr_index+2});

%Find the sequence type
sequence_index=find(ismember(sparheader, 'examination_name')==1);
sequence=strtrim(sparheader{sequence_index+2});



%Philips sdat data has coil elements already combined:
Ncoils=1;

%Now create a record of the dimensions of the data array.  
dims.t=1;
dims.coils=0;
dims.averages=2;
if subspecs>1
    dims.subSpecs=3;
else
    dims.subSpecs=0;
end

specs=fftshift(ifft(fids,[],dims.t),dims.t);

%Now get relevant scan parameters:*****************************

%Leave date blank
date_index=find(ismember(sparheader, 'scan_date')==1);
date = str2num(sparheader{date_index+2});

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
dims.extras=0;



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

%****************************************************************


%Calculate t and ppm arrays using the calculated parameters:
f=[(-spectralwidth/2)+(spectralwidth/(2*sz(1))):spectralwidth/(sz(1)):(spectralwidth/2)-(spectralwidth/(2*sz(1)))];
ppm=-f/(Bo*42.577);
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
out.flags.averaged=0;
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
