%io_loadspec_twix.m
%Jamie Near, McGill University 2014.
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

function out=io_loadspec_twix(filename);


%read in the data using the new mapVBVD.  This code has been adapted to 
%handle both single RAID files and multi-RAID files.  The vast majority of
%Siemens twix data comes as a single RAID file, but I've encoundered a few 
%multi-RAID files, particularly when using VD13D.  The way to distinguish
%them here is that a for a single RAID file, mapVBVD will output a struct, 
%whereas for a multi-RAID file, mapVBVD will output a cell array of structs.
%This code assumes that the data of interest is in the last element of the 
%cell array (possibly a bad assumption under some circumstances):
twix_obj=mapVBVD(filename);
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
sqzDims=twix_obj.image.sqzDims;
sqzSize=twix_obj.image.sqzSize;

%find out if the data was acquired using the rm_special sequence, which 
%does not separate the subspectra into separate array dimensions.  Make 
%changes accordingly.
fid=fopen(filename);
line=fgets(fid);
index=findstr(line,'tSequenceFileName');
equals_index=findstr(line,'= ');
while isempty(index) || isempty(equals_index)
    line=fgets(fid);
    index=findstr(line,'tSequenceFileName');
    equals_index=findstr(line,'= ');
end
sequence=line(equals_index+2:end-1);

%Try to find out what sequnece this is:
isSpecial=~isempty(findstr(sequence,'rm_special')); %Is this Ralf Mekle's SPECIAL Sequence?
isWIP529=~isempty(findstr(sequence,'edit_529')); %Is this WIP 529 (MEGA-PRESS)?
isjnseq=~isempty(findstr(sequence,'jn_')); %Is this one of Jamie Near's sequences?
isMinnMP=~isempty(findstr(sequence,'eja_svs_mpress')); %Is this Eddie Auerbach's MEGA-PRESS?

fclose(fid);
%Ralf Mekle's SPECIAL sequence does not store the subspectra along a
%separate dimension of the data array, so we will separate them
%artifically:
if isSpecial 
    squeezedData=squeeze(dOut.data);
    if twix_obj.image.NCol>1 && twix_obj.image.NCha>1
        data(:,:,:,1)=squeezedData(:,:,[1:2:end-1]);
        data(:,:,:,2)=squeezedData(:,:,[2:2:end]);
        sqzSize=[sqzSize(1) sqzSize(2) sqzSize(3)/2 2];
    elseif twix_obj.NCol>1 && twixObj.image.NCha==1
        data(:,:,1)=squeezedData(:,[1:2:end-1]);
        data(:,:,2)=squeezedData(:,[2:2:end]);
        sqzSize=[sqzSize(1) sqzSize(2)/2 2];
    end
    sqzDims{end+1}='Ida';
else
    data=dOut.data;
end

%Make a pulse sequence identifier for the header (out.seq);
seq=sequence;
  
%Squeeze the data to remove singleton dims
fids=squeeze(data);

%Find the magnetic field strength:
fid=fopen(filename);
line=fgets(fid);
index=findstr(line,'sProtConsistencyInfo.flNominalB0');
equals_index=findstr(line,'= ');
while isempty(index) || isempty(equals_index)
    line=fgets(fid);
    index=findstr(line,'sProtConsistencyInfo.flNominalB0');
    equals_index=findstr(line,'= ');
end
Bo=line(equals_index+1:end);
Bo=str2double(Bo);
fclose(fid);

%Find the number of averages:
fid=fopen(filename);
line=fgets(fid);
index=findstr(line,'<ParamLong."lAverages">');
if isempty(index);
    instancesFound=0;
else
    instancesFound=1;
end
while isempty(index) || instancesFound<RaidLength
    line=fgets(fid);
    index=findstr(line,'<ParamLong."lAverages">');
    if ~isempty(index)
        instancesFound=instancesFound+1;
    end
end
line=fgets(fid);
line=fgets(fid);
Naverages=line;
Naverages=str2num(Naverages);
fclose(fid);

%Find out if multiple coil elements were used:
fid=fopen(filename);
line=fgets(fid);
coil_index=findstr(line,'ParamLong."iMaxNoOfRxChannels"');
while isempty(coil_index)
    line=fgets(fid);
    coil_index=findstr(line,'ParamLong."iMaxNoOfRxChannels"');
end
line=fgets(fid);
line=fgets(fid);
Ncoils=line;
Ncoils=str2num(Ncoils);
fclose(fid);

%Find the TE:
fid=fopen(filename);
line=fgets(fid);
index=findstr(line,'alTE[0]');
equals_index=findstr(line,'= ');
while isempty(index) || isempty(equals_index)
    line=fgets(fid);
    index=findstr(line,'alTE[0]');
    equals_index=findstr(line,'= ');
end
TE=line(equals_index+1:end);
TE=str2double(TE);
fclose(fid);

%Find the TR:
fid=fopen(filename);
line=fgets(fid);
index=findstr(line,'alTR[0]');
equals_index=findstr(line,'= ');
while isempty(index) || isempty(equals_index)
    line=fgets(fid);
    index=findstr(line,'alTR[0]');
    equals_index=findstr(line,'= ');
end
TR=line(equals_index+1:end);
TR=str2double(TR);
fclose(fid);


%Naverages
%Ncoils
%Now begin indexing the dimensions of the data array. ie. create the dims
%structure, which specifies which dimensions of the data array are being
%used to hold the time-domain data, the multiple coil channels, the
%average, the sub-spectra, and any additional dimensions.
sqzDims_update=sqzDims;
dimsToIndex=[1:length(sqzDims)];

%First index the dimension of the time-domain data
dims.t=find(strcmp(sqzDims,'Col'));
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
    dimsToIndex=dimsToIndex(dimsToIndex~=dims.averages);
else
    dims.averages=0;
end

%Now we have indexed the dimensions containing the timepoints, the coil
%channels, and the averages.  As we indexed each dimension, we removed the
%corresponding index from the dimsToIndex vector.  At this point, if there
%are any values left in the dimsToIndex vector, then there must be some
%additional dimensions that need indexing.  We assume that if sub-spectra exist,
%then these must be indexed in either the 'Ida' dimension (for all Jamie
%Near's VB-version pulse sequences), or the 'Eco' dimension (for the WIP
%MEGA-PRESS seuqence or the Minnesota MEGA-PRESS sequence). 
if ~isempty(dimsToIndex)
    %Now index the dimension of the sub-spectra
    if isjnseq  || isSpecial
        if strcmp(version,'vd') || strcmp(version,'ve')
            dims.subSpecs=find(strcmp(sqzDims,'Set'));
        else
            dims.subSpecs=find(strcmp(sqzDims,'Ida'));
        end
    elseif isWIP529 || isMinnMP
        dims.subSpecs=find(strcmp(sqzDims,'Eco'));
    else
        dims.subSpecs=dimsToIndex(1);
    end
    if ~isempty(dims.subSpecs)
        %remove the sub-spectra dimension from the dimsToIndex vector
        dimsToIndex=dimsToIndex(dimsToIndex~=dims.subSpecs);
    else
        dims.subSpecs=0;
    end
else
    dims.subSpecs=0;
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
%   3) averages.
%   4) subSpecs.
%   5) extras.
if length(sqzDims)==5
    fids=permute(fids,[dims.t dims.coils dims.averages dims.subSpecs dims.extras]);
    dims.t=1;dims.coils=2;dims.averages=3;dims.subSpecs=4;dims.extras=5;
elseif length(sqzDims)==4
    if dims.extras==0
        fids=permute(fids,[dims.t dims.coils dims.averages dims.subSpecs]);
        dims.t=1;dims.coils=2;dims.averages=3;dims.subSpecs=4;dims.extras=0;
    elseif dims.subSpecs==0
        fids=permute(fids,[dims.t dims.coils dims.averages dims.extras]);
        dims.t=1;dims.coils=2;dims.averages=3;dims.subSpecs=0;dims.extras=4;
    elseif dims.averages==0
        fids=permute(fids,[dims.t dims.coils dims.subSpecs dims.extras]);
        dims.t=1;dims.coils=2;dims;averages=0;dims.subSpecs=3;dims.extras=4;
    elseif dims.coils==0
        fids=permute(fids,[dims.t dims.averages dims.subSpecs dims.extras]);
        dims.t=1;dims.coils=0;dims.averages=2;dims.subSpecs=3;dims.extras=4;
    end
elseif length(sqzDims)==3
    if dims.extras==0 && dims.subSpecs==0
        fids=permute(fids,[dims.t dims.coils dims.averages]);
        dims.t=1;dims.coils=2;dims.averages=3;dims.subSpecs=0;dims.extras=0;
    elseif dims.extras==0 && dims.averages==0
        fids=permute(fids,[dims.t dims.coils dims.subSpecs]);
        dims.t=1;dims.coils=2;dims.averages=0;dims.subSpecs=3;dims.extras=0;
    elseif dims.extras==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.averages dims.subSpecs]);
        dims.t=1;dims.coils=0;dims.averages=2;dims.subSpecs=3;dims.extras=0;
    end
elseif length(sqzDims)==2
    if dims.extras==0 && dims.subSpecs==0 && dims.averages==0
        fids=permute(fids,[dims.t dims.coils]);
        dims.t=1;dims.coils=2;dims.averages=0;dims.subSpecs=0;dims.extras=0;
    elseif dims.extras==0 && dims.subSpecs==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.averages]);
        dims.t=1;dims.coils=0;dims.averages=2;dims.subSpecs=0;dims.extras=0;
    elseif dims.extras==0 && dims.averages==0 && dims.coils==0
        fids=permute(fids,[dims.t dims.subSpecs]);
        dims.t=1;dims.coils=0;dims.averages=0;dims.subSpecs=2;dims.extras=0;
    end
elseif length(sqzDims)==1
    fids=permute(fids,[dims.t]);
    dims.t=1;dims.coils=0;dims.averages=0;dims.subSpecs=0;dims.extras=0;
end

%Now get the size of the data array:
sz=size(fids);

%Now take fft of time domain to get fid:
specs=fftshift(ifft(fids,[],dims.t),dims.t);
    

%Now get relevant scan parameters:*****************************

%Get Spectral width and Dwell Time
fid=fopen(filename);
line=fgets(fid);
index=findstr(line,'sRXSPEC.alDwellTime[0]');
equals_index=findstr(line,'= ');
if isempty(index) || isempty(equals_index);
    instancesFound=0;
elseif ~isempty(index) && ~isempty(equals_index)
    instancesFound=1;
end
while (isempty(index) || isempty(equals_index)) || instancesFound<RaidLength
    line=fgets(fid);
    index=findstr(line,'sRXSPEC.alDwellTime[0]');
    equals_index=findstr(line,'= ');
    if ~isempty(index) && ~isempty(equals_index)
        instancesFound=instancesFound+1;
    end
end
dwelltime=line(equals_index+1:end);
dwelltime=str2double(dwelltime)*1e-9;
spectralwidth=1/dwelltime;
fclose(fid);
    
%Get TxFrq
fid=fopen(filename);
line=fgets(fid);
index=findstr(line,'sTXSPEC.asNucleusInfo[0].lFrequency');
equals_index=findstr(line,'= ');
while isempty(index) || isempty(equals_index)
    line=fgets(fid);
    index=findstr(line,'sTXSPEC.asNucleusInfo[0].lFrequency');
    equals_index=findstr(line,'= ');
end
txfrq=line(equals_index+1:end);
txfrq=str2double(txfrq);
fclose(fid);

%Get Date
fid=fopen(filename);
line=fgets(fid);
index=findstr(line,'ParamString."atTXCalibDate">');
%quotes_index=findstr(line,'  "');
while isempty(index) %|| isempty(equals_index)
    line=fgets(fid);
    index=findstr(line,'ParamString."atTXCalibDate">');
    if ~isempty(index)
        line=fgets(fid);
        line=fgets(fid);
        quote_index=findstr(line,'  "');
    end
end
date=line(quote_index+3:quote_index+10);
date=str2double(date);
fclose(fid);

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
ppm=-f/(Bo*42.577);
ppm=ppm+4.6082;

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
out.seq=seq;
out.te=TE/1000;
out.tr=TR/1000;
out.pointsToLeftshift=twix_obj.image.freeParam(1);


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



%DONE
