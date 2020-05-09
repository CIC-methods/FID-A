%io_loadspec_twix.m
%Jamie Near, McGill University 2014.
%Edits from Franck Lamberton, 2017.
%
% USAGE:
% out=io_loadspec_twix(filename);
% 
% DESCRIPTION:
% Reads in siemens twix raw data (.dat file) using the mapVBVD.m and 
% twix_map_obj.m functions from Philipp Ehses (philipp.ehses@tuebingen.mpg.de).
% 
% op_loadspec_twix outputs the data in structure format, with fields corresponding to time
% scale, fids, frequency scale, spectra, and header fields containing
% information about the acquisition.  The resulting matlab structure can be
% operated on by the other functions in this MRS toolbox.
% 
% INPUTS:
% filename   = filename of Siemens twix data to load.
%
% OUTPUTS:
% out        = Input dataset in FID-A structure format.

function out=io_loadCSI_twix(filename)


%read in the data using the new mapVBVD.  This code has been adapted to 
%handle both single RAID files and multi-RAID files.  The vast majority of
%Siemens twix data comes as a single RAID file, but I've encoundered a few 
%multi-RAID files, particularly when using VD13D.  The way to distinguish
%them here is that a for a single RAID file, mapVBVD will output a struct, 
%whereas for a multi-RAID file, mapVBVD will ou

%This code assumes that the data of interest is in the last element of the 
%cell array (possibly a bad assumption under some circumstances):
if(~exist('filename', 'var'))
    twix_obj = mapVBVD;
else
    twix_obj=mapVBVD(filename);
end
if isstruct(twix_obj)
    disp('single RAID file detected.');
    RaidLength=1;
elseif iscell(twix_obj)
    disp('multi RAID file detected.');
    RaidLength=length(twix_obj);
    %assume that the data of interest is in the last element of the cell.
    twix_obj=twix_obj{RaidLength};
end
dOut.data=twix_obj.image();
version=twix_obj.image.softwareVersion;
sqzSize=twix_obj.image.sqzSize; 
sqzDims=twix_obj.image.sqzDims;


%find out what sequence, the data were acquired with.  If this is a
%multi-raid file, then the header may contain multiple instances of
%'tSequenceFileName' for different scans (including a pre-scan).
%Therefore, if multi-raid file, we will need to do a bit of extra digging 
%to find the correct sequence name.  
 sequence=twix_obj.hdr.Config.SequenceFileName;  

 isSiemens=(contains(sequence,'svs_se') ||... %Is this the Siemens PRESS seqeunce?
            contains(sequence,'svs_st')) && ... % or the Siemens STEAM sequence?
             ~contains(sequence,'eja_svs');    %And make sure it's not 'eja_svs_steam'.


%Squeeze the data to remove singleton dims
fids=squeeze(dOut.data);

%Make a pulse sequence identifier for the header (out.seq);
seq=sequence;

%Find the magnetic field strength:
Bo=twix_obj.hdr.Dicom.flMagneticFieldStrength;

%Find the number of averages:
Naverages=twix_obj.hdr.Meas.Averages;

%Find out if multiple coil elements were used:
Ncoils=twix_obj.hdr.Meas.iMaxNoOfRxChannels;  

%Find the TE:
TE = twix_obj.hdr.MeasYaps.alTE{1};  %Franck Lamberton

%Find the TR:
TR = twix_obj.hdr.MeasYaps.alTR{1};  %Franck Lamberton

%Now begin indexing the dimensions of the data array. ie. create the dims
%structure, which specifies which dimensions of the data array are being
%used to hold the time-domain data, the multiple coil channels, the
%average, the sub-spectra, and any additional dimensions.
dimsToIndex = 1:length(sqzDims);

%First index the dimension of the time-domain data
dims.t = find(strcmp(sqzDims,'Col'));
if ~isempty(dims.t)
    %remove the time dimension from the dimsToIndex vector
    dimsToIndex=dimsToIndex(dimsToIndex~=dims.t);
else
    dims.t=0;
    error('ERROR:  Spectrom contains no time domain information!!');
end

%Now index the dimension of the coil channels
dims.coils=find(strcmp(sqzDims,'Cha'));
if ~isempty(dims.coils)
    %remove the coils dimension from the dimsToIndex vector
    dimsToIndex=dimsToIndex(dimsToIndex~=dims.coils);
else
    dims.coils=0;
end

%Now index the dimension of the averages
if strcmp(version,'vd') || strcmp(version,'ve')
    dims.averages=find(strcmp(sqzDims,'Ave'));
else
    dims.averages=find(strcmp(sqzDims,'Set'));
end

if ~isempty(dims.averages)
    %remove the averages dimension from the dimsToIndex vector
    dimsToIndex=dimsToIndex(dimsToIndex ~= dims.averages);
else
    %If no Averages dimension was found, then check for a "Repetitions"
    %dimension.  If that is found, store it under "averages".  If both
    %"Averages" and "Repetitions" dimensions are found, "Repetitions" will
    %be indexed under "Extras", since "Repetitions is not currently an
    %option in FID-A.
    dims.averages=find(strcmp(sqzDims,'Rep'));
    if ~isempty(dims.averages)
        dimsToIndex=dimsToIndex(dimsToIndex~=dims.averages);
    else
        %If neither an "Averages" or a "Repetitions" dimension is found,
        %then set the FID-A "Averages" dimension to zero.
        dims.averages=0;
    end
end

%k space coordinates in the x direction 
dims.x = find(strcmp(sqzDims,'Seg'));
if ~isempty(dims.x)
    %remove the coils dimension from the dimsToIndex vector
    dimsToIndex=dimsToIndex(dimsToIndex~=dims.x);
else
    dims.x=0;
end

%k space coordinates in the y direction
dims.y = find(strcmp(sqzDims,'Lin'));
if ~isempty(dims.y)
    %remove the coils dimension from the dimsToIndex vector
    dimsToIndex=dimsToIndex(dimsToIndex~=dims.y);
else
    dims.y=0;
end
%And if any further dimensions exist after indexing the sub-spectra, call
%these the 'extras' dimension.  
if ~isempty(dimsToIndex)
    %Now index the 'extras' dimension
    dims.extras=dimsToIndex(1);
    if ~isempty(dims.extras)
        %remove the extras dimension from the dimsToIndex vector
        dimsToIndex=dimsToIndex(dimsToIndex~=dims.extras);
    else
        dims.extras=0;
    end
else
    dims.extras=0;
end

%Now that we've indexed the dimensions of the data array, we now need to
%permute it so that the order of the dimensions is standardized:  we want
%the order to be as follows:  
%   1) time domain data.  
%   2) coils.
%   4) k_x.
%   5) k_y.
%   6) extras
dimsArray = [dims.t dims.coils, dims.x dims.y dims.extras];
dimsCell = ["t", "coils", "x", "y", "extras"];
for i = 1:numel(dimsArray)
    if(dimsArray(i) == 0)
        dimsArray(i) = [];
        dimsCell(i) = [];
    end
end
fids = permute(fids, dimsArray);
for i = 1:size(dimsArray, 2)
    dims.(dimsCell(i)) = i;
end

%Now get the size of the data array:
sz=size(fids);

%Get Spectral width and Dwell Time
dwelltime = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1}*1e-9;  %Franck Lamberton
spectralwidth=1/dwelltime;

%Get FoV of the CSI image
fovX = twix_obj.hdr.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV;
fovY = twix_obj.hdr.MeasYaps.sSliceArray.asSlice{1}.dPhaseFOV;

%Get TxFrq
txfrq=twix_obj.hdr.Meas.Frequency;

if isSiemens
    leftshift = twix_obj.image.freeParam(1);
else
    leftshift = twix_obj.image.freeParam(1);
end

%****************************************************************

t=0:dwelltime:(sz(1)-1)*dwelltime;

%FILLING IN DATA STRUCTURE
out.fids=fids;
out.sz=sz;
out.t=t;    
out.spectralwidth=spectralwidth;
out.dwelltime=dwelltime;
out.txfrq=txfrq;
out.date=date;
out.dims=dims;
out.Bo=Bo;
out.seq=seq;
out.te=TE/1000;
out.tr=TR/1000;
out.pointsToLeftshift=leftshift;
out.fovX = fovX;    
out.fovY = fovY;    
out.deltaX = fovX/size(fids, dims.x);
out.deltaY = fovY/size(fids, dims.y);
out.deltaK_X = 1/fovX;
out.deltaK_Y = 1/fovY;
out.fovK_X = 1/out.deltaX;
out.fovK_Y = 1/out.deltaY;
out.imageOrigin = [twix_obj.hdr.Config.VoI_Position_Sag, twix_obj.hdr.Config.VoI_Position_Cor, twix_obj.hdr.Config.VoI_Position_Tra];
out.k_XCoordinates = -out.fovK_X/2+out.deltaK_X/2:out.deltaK_X:out.fovK_X/2-out.deltaK_X/2;
out.k_YCoordinates = -out.fovK_Y/2+out.deltaK_Y/2:out.deltaK_Y:out.fovK_Y/2-out.deltaK_Y/2;

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
out.flags.spatialFT = 0;
out.flags.spectralFT = 0;
out.flags.coilCombined = 0;
end

%DONE
