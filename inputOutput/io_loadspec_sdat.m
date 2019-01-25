%io_loadspec_sdat.m
%Jamie Near, McGill University 2014.
%Georg Oeltzschner, Johns Hopkins University 2018.
%
% USAGE:
% out=io_loadspec_sdat(filename,subspecs);
% 
% DESCRIPTION:
% Reads in Philips MRS data (.spar and .sdat files) using code adapted from 
% PhilipsRead.m, provided as part of the Gannet software package by Richard 
% Edden (gabamrs.com).
% 
% io_loadspec_sdat outputs the data in structure format, with fields 
% corresponding to time scale, fids, frequency scale, spectra, and header 
% fields containing information about the acquisition.  The resulting 
% matlab structure can be operated on by the other functions in this MRS 
% toolbox. ALSO:  This code is not currently smart enough to parse out all 
% of the relevant information from the header file, such as the number of 
% subspectra.  So for now, these details must be passed to the function as 
% input arguments.  Help implementing these improvements are most welcome!!
% 
% INPUTS:
% filename   = filename of Philips sdat file to be loaded.
% subspecs   = number of subspectra in the data (from spectral editing, ISIS, etc.)
%
% OUTPUTS:
% out        = Input dataset in FID-A structure format.

function out = io_loadspec_sdat(filename,subspecs)

% Read in the data and header information
[data, header] = philipsLoad(filename);

% Determine the dimensions of the data
data_size       = size(data);
dims.t          = find(data_size == header.samples);
dims.averages   = find(data_size == header.rows);
dims.coils      = 0; % SDAT is already coil-combined
% Now arrange in the standard order (samples-avgs-subspecs):
data = permute(data ,[dims.t dims.averages]);
dims.t = 1;
dims.averages = 2;
dims.extras = 0;

% We have no way of actually knowing the number of sub-spectra (e.g. for 
% MEGA-PRESS or HERMES acquisitions, so we will split the averages 
% according to the 'subspecs' input.
% Initialize fids array:
fids = squeeze(zeros(header.samples, header.rows/subspecs, subspecs));
if subspecs == 2
    %Split the subspectra out of the "averages" dimension:
    fids(:,:,1) = data(:,[1:2:end]);
    fids(:,:,2) = data(:,[2:2:end]);
elseif subspecs == 4
    fids(:,:,1) = data(:,[1:4:end]);
    fids(:,:,2) = data(:,[2:4:end]);
    fids(:,:,3) = data(:,[3:4:end]);
    fids(:,:,4) = data(:,[4:4:end]);
else
    fids = data;
end

if subspecs > 1
    dims.subSpecs = 3;
else
    dims.subSpecs = 0;
end

sz = size(fids);

% Fill in the header information
txfrq = header.synthesizer_frequency; % transmitter frequency [Hz]
Bo = txfrq/42577000; % B0 [T]
averages = size(fids,2)*size(fids,3); % number of rows in the file
rawAverages = averages;
spectralwidth = header.sample_frequency; % bandwidth [Hz]
dwelltime = 1/spectralwidth; % dwelltime [s]
te = header.echo_time; % echo time [ms]
tr = header.repetition_time; % repetition time [ms]
sequence = header.scan_id; % sequence type
Ncoils = 1; % SDAT already coil-combined
rawSubspecs = subspecs;
date = header.scan_date; % scan date [yyyy.mm.dd hh:mm:ss]

% Produce specs
specs = fftshift(ifft(fids,[],dims.t),dims.t);
% Calculate t and ppm arrays using the calculated parameters:
f = [(-spectralwidth/2)+(spectralwidth/(2*sz(1))):spectralwidth/(sz(1)):(spectralwidth/2)-(spectralwidth/(2*sz(1)))];
ppm = -f/(Bo*42.577);
ppm = ppm+4.65;
t = [0:dwelltime:(sz(1)-1)*dwelltime];

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
