%io_CSIload_twix.m
%Brenden Kadota, Jamie Near, Sunnybrook 2021.
%
% USAGE:
% out = io_CSIload_twix('filename');
% out = io_CSIload_twix(scan_number);
%
% DESCRIPTION:
% Reads in siemens twix raw data (.dat file) using the mapVBVD.m and
% twix_map_obj.m functions from Philipp Ehses (philipp.ehses@tuebingen.mpg.de).
%
% io_CSIload_twix takes a CSI scan and outputs a FID-A csi strucure needed
% for later processing steps. Fids dimensiosn are arranged acording to the
% dimension attribute. Image origin in is (x, y, z).
%
% INPUTS:
% filename   = filename or scan number of Siemens twix data to load.
%
% OUTPUTS:
% out        = Input dataset in FID-A structure format.

function MRSIStruct = io_CSIload_twix(filename)
    arguments
        filename (1,:) {mustBeFile}
    end

    %This code assumes that the data of interest is in the last element of the cell array
    twix_obj = readTwixFile(filename);
    dOut.data = twix_obj.image();
    %Squeeze the data to remove singleton dims
    data = squeeze(dOut.data);
    % reverse the direction of rotation
    data = conj(data);

    %get version number and dims
    %version = twix_obj.image.softwareVersion;
    sqzDims = twix_obj.image.sqzDims;

    %find out what sequence, the data were acquired with.
    sequence = twix_obj.hdr.Config.SequenceFileName;

    %Try to find out what sequnece this is:
    %isSiemens = (contains(sequence,'csi_se') ||... %Or the Siemens CSI PRESS sequence?
    %contains(sequence,'csi_st'));

    isRosette = (contains(sequence, 'ros', 'IgnoreCase',true));
    isCartesian = true;
    %%Add more when more non cartesians come along
    if(isRosette)
        isCartesian = false;
    end

    %fill dims field based on dimensions present
    dims = fillDimsField(sqzDims);
    %premute dims to be in a standardized order
    [dims, data] = permuteDims(dims, data);
    %get matrix size
    [numX, numY, numZ] = getNumberOfVoxels(isCartesian, data, dims);


    sz = size(data);
    %Find the number of averages.  'averages' will specify the current numbe
    averages = 1;
    rawAverages = 1;
    if dims.averages ~= 0
        averages = sz(dims.averages);
        rawAverages = averages;
    end

    %Find the number of subspecs.  'subspecs' will specify the current number
    if dims.subSpecs ~= 0
        subspecs = sz(dims.subSpecs);
        rawSubspecs = subspecs;
    else
        subspecs = 1;
        rawSubspecs = subspecs;
    end

    leftshift = twix_obj.image.freeParam(1);
    %Get Spectral width and Dwell Time
    dwelltime = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1}*1e-9;  %Franck Lamberton
    spectralwidth = 1/dwelltime;
    adcTime = 0:dwelltime:(sz(dims.t)-1)*dwelltime;
    if(isCartesian)
        data = flip(data, dims.ky);
        data = flip(data, dims.kx);
    end
    %Now get the size of the data array:
    sz = size(data);
    %****************************************************************
    %FILLING IN DATA STRUCTURE
    MRSIStruct.data = data;
    MRSIStruct.sz = sz;
    MRSIStruct.adcTime = adcTime;
    if(isCartesian)
        MRSIStruct.spectralWidth = spectralwidth; 
        MRSIStruct.spectralTime = adcTime;
        MRSIStruct.spectralDwellTime = dwelltime;
    end
    MRSIStruct.adcDwellTime = dwelltime;
    MRSIStruct.txfrq = twix_obj.hdr.Meas.lFrequency;
    MRSIStruct.scanDate = findScanDate(twix_obj);
    MRSIStruct.dims = dims;
    MRSIStruct.Bo = twix_obj.hdr.Dicom.flMagneticFieldStrength;
    MRSIStruct.seq = sequence;
    MRSIStruct.te = twix_obj.hdr.MeasYaps.alTE{1}/1000;
    MRSIStruct.tr = twix_obj.hdr.MeasYaps.alTR{1}/1000;
    MRSIStruct.pointsToLeftshift = twix_obj.image.freeParam(1);
    MRSIStruct = findAndSetFov(MRSIStruct, twix_obj);
    MRSIStruct = calcualteVoxelSize(MRSIStruct, numX, numY, numZ);
    MRSIStruct.averages = averages;
    MRSIStruct.rawAverages = rawAverages;
    MRSIStruct.subspecs = subspecs;
    MRSIStruct.rawSubspecs = rawSubspecs;
    MRSIStruct.pointsToLeftshift = leftshift;
    MRSIStruct = findImageOrigin(MRSIStruct, twix_obj);
    MRSIStruct = calculateVoxelCoodinates(MRSIStruct);
    MRSIStruct = calculateAffineMatrix(MRSIStruct, twix_obj);
    MRSIStruct = setDefaultFlagValues(MRSIStruct, isCartesian);
end


%convert rotation around vector to rotation matrix. Formula found from
%wikipedia:
%https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle.
function rotation_matrix = getRotationMatrixFromVector(x, y, z, theta)
    vect = [x,y,z];
    rotation_matrix = zeros(3,3);
    rotation_matrix(1,1) = cos(theta)+vect(1)^2*(1-cos(theta));
    rotation_matrix(1,2) = vect(1)*vect(2)*(1-cos(theta))-vect(3)*sin(theta);
    rotation_matrix(1,3) = vect(1)*vect(3)*(1-cos(theta))+vect(2)*sin(theta);
    rotation_matrix(2,1) = vect(2)*vect(1)*(1-cos(theta))+vect(3)*sin(theta);
    rotation_matrix(2,2) = cos(theta)+vect(2)^2*(1-cos(theta));
    rotation_matrix(2,3) = vect(2)*vect(3)*(1-cos(theta))-vect(1)*sin(theta);
    rotation_matrix(3,1) = vect(3)*vect(1)*(1-cos(theta))-vect(2)*sin(theta);
    rotation_matrix(3,2) = vect(3)*vect(2)*(1-cos(theta))+vect(2)*sin(theta);
    rotation_matrix(3,3) = cos(theta) + vect(3)^2*(1-cos(theta));
end

function out = initalizeZeroIfEmpty(value)
    if(isempty(value))
        out = 0;
    else
        out = value;
    end
end

function out = initalizeOneIfEmpty(value)
    if(isempty(value))
        out = 1;
    else
        out = value;
    end
end

function out = setDefaultFlagValues(out, isCartesian)
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

function twix_obj = readTwixFile(filename)
    if(~exist('filename', 'var'))
        twix_obj = mapVBVD;
    else
        twix_obj = mapVBVD(char(filename));
    end

    if isstruct(twix_obj)
        disp('single RAID file detected.');
    elseif iscell(twix_obj)
        disp('multi RAID file detected.');
        RaidLength = length(twix_obj);
        %assume that the data of interest is in the last element of the cell.
        twix_obj = twix_obj{RaidLength};
    end
end

function [z_vect, theta] = getZVectorAndTheta(twix_obj)
    %initalize to zero if empty
    if(isfield(twix_obj.hdr.Meas, 'Meas.VoI_Normal_Sag'))
        z_vect(1) = -initalizeZeroIfEmpty(twix_obj.hdr.Meas.VoI_Normal_Sag);
        z_vect(2) = -initalizeZeroIfEmpty(twix_obj.hdr.Meas.VoI_Normal_Cor);
        z_vect(3) = initalizeOneIfEmpty(twix_obj.hdr.Meas.VoI_Normal_Tra);
        theta = initalizeZeroIfEmpty(twix_obj.hdr.Meas.VoI_InPlaneRotAngle);
    else
        z_vect(1) = -initalizeZeroIfEmpty(twix_obj.hdr.Meas.VoiNormalSag);
        z_vect(2) = -initalizeZeroIfEmpty(twix_obj.hdr.Meas.VoiNormalCor);
        z_vect(3) = initalizeOneIfEmpty(twix_obj.hdr.Meas.VoiNormalTra);
        theta = initalizeZeroIfEmpty(twix_obj.hdr.Meas.VoiInPlaneRot);
    end
end

function MRSIStruct = calculateAffineMatrix(MRSIStruct, twixObj)
    [zVector, theta] = getZVectorAndTheta(twixObj);
    rotationMatrix = getRotationMatrixFromVector(zVector(1), zVector(2), zVector(3), theta);
    %create affine matrix from rotation matrix
    rotationMatrix(4,4) = 1;

    %find y vector perpendicular to z with x equal to zero
    yVector = [0, 1, -zVector(2)/zVector(3)];
    %rotate y vector around z using rotation matrix
    yVector = yVector/norm(yVector);
    yVector = rotationMatrix*[yVector'; 1];
    yVector = yVector(1:3);
    %get x vector from taking the cross product of y and z. Ensures x y z are
    %perpendicular
    xVector = cross(yVector, zVector);

    %get affine rotation matrix
    affineRotationMatrix = [xVector', yVector, zVector', [0,0,0]'; 0 0 0 1];

    %get scaling matrix
    affineScaleMatrix = eye(4);
    affineScaleMatrix(1,1) = MRSIStruct.voxelSize.x;
    affineScaleMatrix(2,2) = MRSIStruct.voxelSize.y;
    affineScaleMatrix(3,3) = MRSIStruct.voxelSize.z;

    %get translate matrix
    affineTranslateMatrix = eye(4);
    affineTranslateMatrix(1,4) = -MRSIStruct.imageOrigin(1) - MRSIStruct.fov.x/2;
    affineTranslateMatrix(2,4) = -MRSIStruct.imageOrigin(2) - MRSIStruct.fov.y/2;
    affineTranslateMatrix(3,4) = MRSIStruct.imageOrigin(3) - MRSIStruct.fov.z/2;

    %ORDER MATTERS HERE!!! First we scale image coordinates to be voxel sizes.
    %Then we shift the voxels so the first voxel starts at smallest x,y,z
    %coordinate. Then we rotate the voxels using the rotation matrix
    affineMatrix = affineRotationMatrix*affineTranslateMatrix*affineScaleMatrix;
    MRSIStruct.affineMatrix = affineMatrix;
end

function MRSIStruct = findImageOrigin(MRSIStruct, twix_obj)
    MRSIStruct.imageOrigin = zeros(1,3);
    if(isfield(twix_obj.hdr.Config, 'VoI_Position_Sag'))
        fields = ["VoI_Position_Sag", "VoI_Position_Cor", "VoI_Position_Tra"];
    elseif(isfield(twix_obj.hdr.Config, 'Voi_Position_Sag'))
        fields = ["Voi_Position_Sag", "Voi_Position_Cor", "Voi_Position_Tra"];
    else
        fields = [];
    end

    for i = 1:length(fields)
        if(~isempty(twix_obj.hdr.Config.(fields(i))))
            MRSIStruct.imageOrigin(i) = twix_obj.hdr.Config.(fields(i));
        end
    end
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

%calculate the number of x and y voxels in the image dimension
function [numX, numY, numZ] = getNumberOfVoxels(isCartesian, data, dims)
    if(~isCartesian)
        numX = -1;
        numY = -1;
        numZ = -1;
        while numX < 0 && numY < 0
            numX = input("Please enter an integer for Nx: \n");
            numY = input("Please enter an integer for Ny: \n");
            numZ = input("Please enter an integer for Nz: \n");
            if(~isa(numX, 'double') || ~isa(numY,'double'))
                disp('Please enter a number')
                continue
            end
            if(numX < 0 || floor(numX) ~= numX || numY < 0 || floor(numY) ~= numY)
                disp('Please enter an integer > 0')
            end
        end
    else
        numX = size(data, dims.kx);
        numY = size(data, dims.ky);
        numZ = 1;
        if(dims.kz ~= 0); numZ = size(data, dims.kz); end
    end
end

function [MRSIStruct] = findAndSetFov(MRSIStruct, twix_obj)
    %Get FoV of the CSI image
    fovX = twix_obj.hdr.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV;
    fovY = twix_obj.hdr.MeasYaps.sSliceArray.asSlice{1}.dPhaseFOV;
    fovZ = twix_obj.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness;
    
    MRSIStruct.fov.x = fovX;
    MRSIStruct.fov.y = fovY;
    MRSIStruct.fov.z = fovZ;
end

%get scan date from header fields
function scanDate = findScanDate(twix_obj)
    %use a regular expression to extract date data
    scanDate = regexp(twix_obj.hdr.MeasYaps.tReferenceImage0,...
        '\.(?<year>\d{4})(?<month>\d{2})(?<day>\d{2})', 'names');
    %set date using datetime type
    scanDate = datetime(str2double(scanDate.year), str2double(scanDate.month), str2double(scanDate.day));
end

%calculate and set voxel size
function MRSIStruct = calcualteVoxelSize(MRSIStruct, numX, numY, numZ)
    MRSIStruct = setVoxelSize(MRSIStruct, 'x', getFov(MRSIStruct, 'x')/numX);
    MRSIStruct = setVoxelSize(MRSIStruct, 'y', getFov(MRSIStruct, 'y')/numY);
    MRSIStruct = setVoxelSize(MRSIStruct, 'z', getFov(MRSIStruct, 'z')/numZ);
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