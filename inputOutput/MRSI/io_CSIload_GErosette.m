% io_CSIload_GErosette.m
% Jamie Near, Sunnybrook 2023.
% Louis Lauzon, University of Calgary 2023.
%
% USAGE:
% out = io_CSIload_GErosette('filename');
%
% DESCRIPTION:
% Reads in a GE Rosette MRSI dataset using the GELoad.m function.  
% 
% As input, this script accepts both .h5 files (for with GE systems with 
% Version MR30), or P-files (for earlier system version). 
%
% io_CSIload_GE takes a rosette MRSI scan and outputs a FID-A csi strucure needed
% for later processing steps. Fids dimensions are arranged acording to the
% dimension attribute. Image origin in is (x, y, z).
%
% INPUTS:
% filename   = Name (and optionally the path) of the p-file, passed as a char
%
% OUTPUTS:
% out        = Input dataset in FID-A structure format.

function out = io_CSIload_GErosette(filename)
    arguments
        filename (1,:) {mustBeFile}
    end

    %This code assumes that the data of interest is in the last element of the cell array
    [data,GEhdr] = GELoad_CSI(filename);

    %record the initial size of the data array:
    sz_init=size(data);

    isRosette = (contains(GEhdr.pat_nam, 'rose', 'IgnoreCase',true) || contains(GEhdr.psd_nam, 'eCSI_rst','IgnoreCase',true));
    isCartesian = true;
    %%Add more when more non cartesians come along
    if(isRosette)
        isCartesian = false;
    end

    %fill dims field based on dimensions present
    dims = fillDimsField(GEhdr.dims);
    %premute dims to be in a standardized order
    [dims, data] = permuteDims(dims, data);
    %get matrix size
    numX=GEhdr.n_x;
    numY=GEhdr.n_x;
    numZ=GEhdr.n_slc;

    sz = size(data);
    
    %Find the number of averages.  'averages' will specify the current numbe
    averages = GEhdr.n_avg;
    rawAverages = averages;

    %Squeeze the averages dimension if Averages is equal to 1:
    if averages==1
        data=squeeze(data);
        if dims.t>dims.averages
            dims.t=dims.t-1;
        end
        if dims.coils>dims.averages;
            dims.coils=dims.coils-1;
        end
        if dims.kx>dims.averages;
            dims.kx=dims.kx-1;
        end
        if dims.ky>dims.averages;
            dims.ky=dims.ky-1;
        end
        if dims.kx>dims.averages;
            dims.kz=dims.kz-1;
        end
        if dims.x>dims.averages;
            dims.x=dims.x-1;
        end
        if dims.y>dims.averages;
            dims.y=dims.y-1;
        end
        if dims.z>dims.averages;
            dims.z=dims.z-1;
        end
        if dims.extras>dims.averages;
            dims.extras=dims.extras-1;
        end
        if dims.subSpecs>dims.averages;
            dims.subSpecs=dims.subSpecs-1;
        end
        dims.averages=0;
    end
    

    %Find the number of subspecs.  'subspecs' will specify the current number
    if dims.subSpecs ~= 0
        subspecs = sz(dims.subSpecs);
        rawSubspecs = subspecs;
    else
        subspecs = 1;
        rawSubspecs = subspecs;
    end

    %Set leftshift to zero for now:
    leftshift = 0;

    %Get Spectral width and Dwell Time
    dwelltime = GEhdr.ptl_ns * GEhdr.ptl_dt * 1e-6;
    spectralwidth = 1/dwelltime;
    adcTime = 0:dwelltime:(sz(dims.t)-1)*dwelltime;
    if(isCartesian)
        data = flip(data, dims.ky);
        data = flip(data, dims.kx);
    end
    %Now get the size of the data array:
    sz = size(data);
    
    %Getting Nucleus - PT,2022
%     nucleus=twix_obj.hdr.Config.Nucleus;
      nucleus='1H';
%     switch nucleus
%     case '1H'
        gamma=42.576;
%     case '31P'
%         gamma=17.235;
%     case '13C'
%         gamma=10.7084;
%     end

    %****************************************************************
    %FILLING IN DATA STRUCTURE
    out.data = data;
    out.sz = sz;
    out.adcTime = adcTime;
    if(isCartesian)
        out.spectralWidth = spectralwidth; 
        out.spectralTime = adcTime;
        out.spectralDwellTime = dwelltime;
    end
    out.adcDwellTime = dwelltime;
    out.txfrq = GEhdr.Bo * gamma * 1e6;
    out.scanDate = GEhdr.scn_date;
    out.dims = dims;
    out.Bo = GEhdr.Bo;
    out.nucleus=nucleus;
    out.gamma=gamma;
    out.seq = GEhdr.ser_dsc;
    out.te = GEhdr.t_ech;
    out.tr = GEhdr.t_rep;
    out.fov.x = GEhdr.fov*10;
    out.fov.y = GEhdr.fov*10;
    out.fov.z = GEhdr.slthk;
    out.voxelSize.x = out.fov.x / GEhdr.n_x;
    out.voxelSize.y = out.fov.y / GEhdr.n_x;
    out.voxelSize.z = out.fov.z / 1;
    out.averages = averages;
    out.rawAverages = rawAverages;
    out.subspecs = subspecs;
    out.rawSubspecs = rawSubspecs;
    out.pointsToLeftshift = leftshift;
    out.imageOrigin(1) = 0; %Set to zero for now;
    out.imageOrigin(2) = 0; %Set to zero for now;
    out.imageOrigin(3) = 0; %Set to zero for now;
    out = calculateVoxelCoodinates(out);
    out.affineMatrix = []; %Leave empty for now;
    
    %FILLING IN THE FLAGS
    out.flags.writtentostruct = 1;
    out.flags.gotparams = 1;
    out.flags.leftshifted = 0;
    out.flags.filtered = 0;
    out.flags.zeropadded = 0;
    out.flags.freqcorrected = 0;
    out.flags.phasecorrected = 0;
    out.flags.averaged = 0;
    out.flags.addedrcvrs = 0;
    out.flags.subtracted = 0;
    out.flags.writtentotext = 0;
    out.flags.downsampled = 0;
    out.flags.spatialFT = 0;
    out.flags.spectralFT = 0;
    out.flags.coilCombined = 0;
    out.flags.isFourSteps = 0;
    out.flags.isCartesian = isCartesian;
end

function dims = fillDimsField(sqzDims)
    siemensLabels = {'Col', 'Cha', 'Ave', 'Set', 'Rep', 'Seg', 'Phs', 'Lin', 'Sli'};
    dimLabelMRSI = {'t', 'coils', 'averages', 'averages', 'averages', 'kx', 'kx', 'ky', 'kz'};
    SiemenstoMRSI = containers.Map(siemensLabels, dimLabelMRSI);
    dims.t = 0;
    dims.coils = 0;
    dims.averages = 0;
    dims.kx = 0;
    dims.ky = 0;
    dims.kz = 0;
    dims.x = 0;
    dims.y = 0;
    dims.z = 0;
    dims.extras = 0;
    dims.subSpecs = 0;
    for i = 1:length(sqzDims)
        if(isKey(SiemenstoMRSI, sqzDims{i}))
            dimLabel = SiemenstoMRSI(sqzDims{i});
            dims.(dimLabel) = i;
        else
            if(dims.extras == 0)
                dims.extras = i;
            else
                error('Two extra fields! Not sure what to do')
            end
        end
    end
end

function [dims, data] = permuteDims(dims, data)
    %Now that we've indexed the dimensions of the data array, we now need to
    %permute it so that the order of the dimensions is standardized:  we want
    %the order to be as follows:
    %   1) time domain data.
    %   2) coils.
    %   3) averages.
    %   4) subSpecs.
    %   5) csi_x.
    %   6) csi_y.
    %   7) extras.
    dimsArray = [dims.t dims.coils, dims.averages dims.subSpecs dims.kx dims.ky dims.kz dims.extras];
    fields = {"t", "coils", "averages", "subSpecs", "kx", "ky", 'kz', "extras"};
    nonZeroIndex = dimsArray ~= 0;
    data = permute(data, dimsArray(nonZeroIndex));

    nonZeroFields = fields(nonZeroIndex);
    for i = 1:length(nonZeroFields)
        dims.(nonZeroFields{i}) = i;
    end
end

%calculate the coordinates of voxels in the image domain
function MRSIStruct = calculateVoxelCoodinates(MRSIStruct)
    %calculate x coordinates
    fovX = getFov(MRSIStruct, 'x');
    voxSizeX = getVoxSize(MRSIStruct, 'x');
    xCoordinates = createCoordinates(fovX/2, voxSizeX);
    xCoordinates = xCoordinates - getImageOrigin(MRSIStruct, 'x');

    %calculate y coordinates
    fovY = getFov(MRSIStruct, 'y');
    voxSizeY = getVoxSize(MRSIStruct, 'y');
    yCoordinates = createCoordinates(fovY/2, voxSizeY);
    yCoordinates = yCoordinates - getImageOrigin(MRSIStruct, 'y');

    %set Coordinates to MRSIStruct
    MRSIStruct = setCoordinates(MRSIStruct, 'x', xCoordinates);
    MRSIStruct = setCoordinates(MRSIStruct, 'y', yCoordinates);
end

