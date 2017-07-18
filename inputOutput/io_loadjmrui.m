%io_loadjmrui.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% out=io_loadjmrui(filename);
% 
% DESCRIPTION:
% Load a jMRUI text file into matlab structure format.  
% 
% INPUTS:
% filename   = filename of the jMRUI txt file.
%
% OUTPUTS:
% out = Input dataset in FID-A structure format.

function out=io_loadjmrui(filename);

%LOAD IN JMRUI .txt FILE
[RF,info]=io_readjmrui(filename);


txfrq=str2num(info.TransmitterFrequency);
sz=[str2num(info.PointsInDataset) 1];
dwelltime=str2num(info.SamplingInterval)*1e-3;
spectralwidth=1/dwelltime;
date=str2num(info.DateOfExperiment);
Bo=txfrq/42577000;  %JMRUI HEADER INCORRECTLY SAYS 3.0

t=[0:dwelltime:sz(1)*dwelltime-dwelltime];

dims.t=1;
dims.coils=0;
dims.averages=0;
dims.subSpecs=0;
dims.extras=0;


fids=RF(:,1);

specs=fftshift(ifft(fids,[],dims.t),dims.t);

f=[(-spectralwidth/2)+(spectralwidth/(2*sz(1))):spectralwidth/(sz(1)):(spectralwidth/2)-(spectralwidth/(2*sz(1)))];
ppm=-f/(Bo*42.577);
ppm=ppm+4.65;

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
out.seq='';
out.te=[];
out.tr=[];
out.pointsToLeftshift=0;

%Write flags
out.flags.writtentostruct=1;
out.flags.gotparams=1;
out.flags.leftshifted=0;
out.flags.filtered=0;
out.flags.zeropadded=0;
out.flags.freqcorrected=0;
out.flags.phasecorrected=0;
out.flags.averaged=1;
out.flags.addedrcvrs=1;
out.flags.writtentotext=0;
out.flags.downsampled=0;
out.flags.isISIS=0;
out.flags.subtracted=1;


