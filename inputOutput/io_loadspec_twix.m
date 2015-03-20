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
    dOut.data=twix_obj.image();
elseif iscell(twix_obj)
    disp('multi RAID file detected.');
    RaidLength=length(twix_obj);
    %assume that the data of interest is in the last element of the cell.
    dOut.data=twix_obj{RaidLength}.image();
end

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
sequence=line(equals_index+1:end);
isSpecial=~isempty(findstr(sequence,'rm_special'));
isWIP529=~isempty(findstr(sequence,'edit_529')); %Is this WIP 529 (MEGA-PRESS)?
fclose(fid);
if isSpecial 
    squeezedData=squeeze(dOut.data);
    data(:,:,:,1)=squeezedData(:,:,[1:2:end-1]);
    data(:,:,:,2)=squeezedData(:,:,[2:2:end]);
else
    data=dOut.data;
end

fids=squeeze(data);

%If we are dealing with WIP 529, then we need to swap the order of the 
%averages and subspecs dimensions:
if isWIP529
    fids=permute(fids,[1,2,4,3]);
end

sz=size(fids);

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


%Naverages
%Ncoils
%Now create a record of the dimensions of the data array.  

if ndims(fids)==5
    dims.t=1;
    dims.coils=2;
    dims.averages=3;
    dims.subSpecs=4;
    dims.extra=5;
elseif ndims(fids)==4  %Default config when 4 dims are acquired
    dims.t=1;
    dims.coils=2;
    dims.averages=3;
    dims.subSpecs=4;
    dims.extra=0;
elseif ndims(fids)<4  %To many permutations...ask user for dims.
    if Naverages == 1 && Ncoils == 1
        if ndims(fids)>1
            dims.t=1;
            dims.coils=0;
            dims.averages=0;
            dims.subSpecs=2;
            dims.extra=0;
        else
            dims.t=1;
            dims.coils=0;
            dims.averages=0;
            dims.subSpecs=0;
            dims.extra=0;
        end
    elseif Naverages>1 && Ncoils==1
        if ndims(fids)>2
            dims.t=1;
            dims.coils=0;
            dims.averages=2;
            dims.subSpecs=3;
            dims.extra=0;
        else
            dims.t=1;
            dims.coils=0;
            dims.averages=2;
            dims.subSpecs=0;
            dims.extra=0;
        end
    elseif Naverages==1 && Ncoils>1
        if ndims(fids)>2
            dims.t=1;
            dims.coils=2;
            dims.averages=0;
            dims.subSpecs=3;
            dims.extra=0;
        else
            dims.t=1;
            dims.coils=2;
            dims.averages=0;
            dims.subSpecs=0;
            dims.extra=0;
        end
    elseif Naverages>1 && Ncoils>1
        if ndims(fids)>3
            dims.t=1;
            dims.coils=2;
            dims.averages=3;
            dims.subSpecs=4;
            dims.extra=0;
        else
            dims.t=1;
            dims.coils=2;
            dims.averages=3;
            dims.subSpecs=0;
            dims.extra=0;
        end
    end
%     dims.t=1;
%     dims.coils=input('Enter the coils Dimension (0 for none):  ');
%     dims.averages=input('Enter the averages Dimension (0 for none):  ');
%     dims.subSpecs=input('Enter the subSpecs Dimension (0 for none);  ');
end

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
if RaidLength==1
    out.pointsToLeftshift=twix_obj.image.freeParam(1);
elseif RaidLength>1
    out.pointsToLeftshift=twix_obj{RaidLength}.image.freeParam(1);
end


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
out.flags.avgNormalized=0;
if out.dims.subSpecs==0
    out.flags.isISIS=0;
else
    out.flags.isISIS=(out.sz(out.dims.subSpecs)==4);
end



%DONE
