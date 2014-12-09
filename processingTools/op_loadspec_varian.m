%op_loadspec_varian.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% out=op_loadspec_varian(filename);
% 
% DESCRIPTION:
% Reads in varian .fid data using the readfid.m and readprocpar.m functions 
% from Martyn Klassen (mklassen@robarts.ca).
% 
% op_loadspec_varian outputs the data in structure format, with fields corresponding to time
% scale, fids, frequency scale, spectra, and header fields containing
% information about the acquisition.  The resulting matlab structure can be
% operated on by the other functions in this MRS toolbox.
% 
% INPUTS:
% filename   = filename of Varian .fid data to load.

function out=op_loadspec_varian(filename);

%read in the data using read_meas_dat
%dOut=op_read_meas_dat2(filename);
%[par,img,k,fids]=fidread(filename);
[fids,hdr,block_hdr]=readfid(filename);
par=readprocpar(filename);


fids=squeeze(fids);
sz=size(fids);

%un-interleave the data for megapress2.c
tempfids(:,[1:1:sz(2)/2],1)=fids(:,[1:2:sz(2)]);
tempfids(:,[1:1:sz(2)/2],2)=fids(:,[2:2:sz(2)]);
fids=squeeze(tempfids);
%fids=fids([end:-1:1],:,:);


sz=size(fids);


%Now create a record of the dimensions of the data array.  
dims.t=1;
dims.coils=0;
dims.averages=2;
dims.subSpecs=3;

specs=fftshift(ifft(fids,[],dims.t),dims.t);

%Now get relevant scan parameters:*****************************

%Get Spectral width and Dwell Time
spectralwidth=par.sw;
dwelltime=1/spectralwidth;

    
%Get TxFrq
txfrq=par.sfrq*1e6;


%Get Date
date=par.date;


%GET THE NUMBER OF RAW AVERAGES AND SUBSPECTRA:

%Find the number of averages.  'averages' will specify the current number
%of averages in the dataset as it is processed, which may be subject to
%change.  'rawAverages' will specify the original number of acquired 
%averages in the dataset, which is unchangeable.
if dims.subSpecs ~=0
    averages=sz(dims.averages)*sz(dims.subSpecs);
    rawAverages=averages;
else
    averages=sz(dims.averages);
    rawAverages=averages;
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
ppm=f/(txfrq/1e6);
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
out.Bo=out.txfrq/42.577/1e6;

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
out.flags.isISIS=0;



%DONE
