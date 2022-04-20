%io_CSIload_rda.m
%Brenden Kadota, SunnyBrook Hosptial 2021.
%
% USAGE:
% [out] = io_CSIload_rda(rda_filename);
%
% DESCRIPTION:
% Reads in siemens rda data (.rda file).
%
% io_CSIload_rda outputs the data in CSI structure format, with fields corresponding to time
% scale, fids, frequency scale, spectra, and header fields containing
% information about the acquisition.  The resulting matlab structure can be
% operated on by the other functions in this CSI toolbox.
%
% INPUTS:
% rda_filename   = filename of Siemens rda data to load.
%
% OUTPUTS:
% out        = Input dataset in FID-A structure format.

function [out] = io_CSIload_rda(rda_filename)

    file = fopen(rda_filename);

    head_start_text = '>>> Begin of header <<<';
    head_end_text   = '>>> End of header <<<';

    tline = fgets(file);

    while (~contains(tline , head_end_text))
        tline = fgets(file);
        if (~contains(tline , head_start_text) + ~contains(tline , head_end_text ) == 2)

            % Store this data in the appropriate format

            occurence_of_colon = strfind(tline, ':');
            variable = tline(1:occurence_of_colon-1) ;
            value    = tline(occurence_of_colon+1 : length(tline)) ;

            switch variable
                case { 'PatientID' , 'PatientName' , 'StudyDescription' , 'PatientBirthDate' , 'StudyDate' , 'StudyTime' , 'PatientAge' , 'SeriesDate' , ...
                        'SeriesTime' , 'SeriesDescription' , 'ProtocolName' , 'PatientPosition' , 'ModelName' , 'StationName' , 'InstitutionName' , ...
                        'DeviceSerialNumber', 'InstanceDate' , 'InstanceTime' , 'InstanceComments' , 'SequenceName' , 'SequenceDescription' , 'Nucleus' ,...
                        'TransmitCoil' }
                    rda.(variable) = value;
                case { 'PatientSex' }
                    % Sex converter! (int to M,F,U)
                    switch value
                        case 0
                            rda.sex = 'Unknown';
                        case 1
                            rda.sex = 'Male';
                        case 2

                            rda.sex = 'Female';
                    end

                case {  'SeriesNumber' , 'InstanceNumber' , 'AcquisitionNumber' , 'NumOfPhaseEncodingSteps' , 'NumberOfRows' , 'NumberOfColumns' , 'VectorSize' }
                    %Integers
                    rda.(variable) = str2double(value);
                case { 'PatientWeight' , 'TR' , 'TE' , 'TM' , 'DwellTime' , 'NumberOfAverages' , 'MRFrequency' , 'MagneticFieldStrength' , 'FlipAngle' , ...
                        'SliceThickness' ,  'FoVHeight' , 'FoVWidth', 'FoV3D' , 'PercentOfRectFoV' , 'PixelSpacingRow' , 'PixelSpacingCol', 'PixelSpacing3D'}
                    %Floats
                    rda.(variable) = str2double(value);
                case {'SoftwareVersion[0]' }
                    rda.software_version = value;
                case {'CSIMatrixSize[0]' }
                    rda.CSIMatrix_Size(1) = str2double(value);
                case {'CSIMatrixSize[1]' }
                    rda.CSIMatrix_Size(2) = str2double(value);
                case {'CSIMatrixSize[2]' }
                    rda.CSIMatrix_Size(3) = str2double(value);
                case {'CSIMatrixSizeOfScan[0]'}
                    rda.CSIMatrixSizeOfScan(1) = str2double(value);
                case {'CSIMatrixSizeOfScan[1]'}
                    rda.CSIMatrixSizeOfScan(2) = str2double(value);
                case {'CSIMatrixSizeOfScan[2]'}
                    rda.CSIMatrixSizeOfScan(3) = str2double(value);
                case {'PositionVector[0]' }
                    rda.PositionVector(1) = str2double(value);
                case {'PositionVector[1]' }
                    rda.PositionVector(2) = str2double(value);
                case {'PositionVector[2]' }
                    rda.PositionVector(3) = str2double(value);
                case {'RowVector[0]' }
                    rda.RowVector(1) = str2double(value);
                case {'RowVector[1]' }
                    rda.RowVector(2) = str2double(value);
                case {'RowVector[2]' }
                    rda.RowVector(3) = str2double(value);
                case {'ColumnVector[0]' }
                    rda.ColumnVector(1) = str2double(value);
                case {'ColumnVector[1]' }
                    rda.ColumnVector(2) = str2double(value);
                case {'ColumnVector[2]' }
                    rda.ColumnVector(3) = str2double(value);
                case {'VOIPositionSag' }
                    rda.imageOrigin(1) = str2double(value);
                case {'VOIPositionCor' }
                    rda.imageOrigin(2) = str2double(value);
                case {'VOIPositionTra' }
                    rda.imageOrigin(3) = str2double(value);
                case {'VOINormalSag' }
                    rda.sliceVector(1) = str2double(value);
                case {'VOINormalCor' }
                    rda.sliceVector(2) = str2double(value);
                case {'VOINormalTra' }
                    rda.sliceVector(3) = str2double(value);

                otherwise
                    % We don't know what this variable is.  Report this just to keep things clear
                    %disp(['Unrecognised variable ' , variable ]);
            end

        else
            % Don't bother storing this bit of the output
        end

    end

    %
    % So now we should have got to the point after the header text
    %
    % Siemens documentation suggests that the data should be in a double complex format (8bytes for real, and 8 for imaginary?)
    %

    bytes_per_point = 16;
    complex_data = fread(file , rda.CSIMatrix_Size(1) * rda.CSIMatrix_Size(2) * rda.CSIMatrix_Size(3) *rda.VectorSize * 2 , 'double');
    %fread(fid , 1, 'double');  %This was a check to confirm that we had read all the data (it passed!)
    fclose(file);

    % Now convert this data into something meaningful

    %Reshape so that we can get the real and imaginary separated
    hmm = reshape(complex_data,  2 , rda.VectorSize , rda.CSIMatrix_Size(1) ,  rda.CSIMatrix_Size(2) ,  rda.CSIMatrix_Size(3) );

    %Combine the real and imaginary into the complex matrix
    fids = complex(hmm(1,:,:,:,:),hmm(2,:,:,:,:));
    fids = squeeze(fids);
    
    %fids = conj(fids);

    % make calculations for the output mrs structure
    sz = size(fids);
    dwelltime = rda.DwellTime/1000000;
    spectralwidth = 1/dwelltime;
    txfrq = rda.MRFrequency*1000000;
    dims.t = 1;
    dims.subSpecs = 0;
    dims.x = 2;
    fids = flip(fids, dims.x);
    dims.y = 3;
    dims.z = 0;
    dims.coils = 0;
    dims.averages = 0;
    dims.extras = 0;
    dims.ky = 0;
    dims.kx = 0;
    Bo = rda.MagneticFieldStrength;
    rawAverages = rda.NumberOfAverages;
    averages = 1;
    subspecs =1;
    rawSubspecs = 'na';
    date = strip(rda.StudyDate);
    seq = rda.SequenceDescription;
    TE = rda.TE;
    TR = rda.TR;
    pointsToLeftShift = 'N/A';

    %Calculate t and ppm arrays using the calculated parameters:
    f = createCoordinates(spectralwidth/2, spectralwidth/sz(1));
    ppm = -f/(Bo*42.577);
    ppm = ppm+4.6082;

    t = 0:dwelltime:(sz(1)-1)*dwelltime;

    out = struct();
    %FILLING IN DATA STRUCTURE
    out = setData(out, fids);
    out = setSize(out, sz);
    out = setPPM(out, ppm);
    out = setAdcTime(out, t);
    out = setSpectralTime(out, t);
    out = setSpectralDwellTime(out, dwelltime);
    out = setSpectralWidth(out, spectralwidth);
    out.txfrq = txfrq;
    out.date = datetime(str2double(date(1:4)), str2double(date(5:6)), str2double(date(7:8)));
    out.dims = dims;
    out.Bo = Bo;
    out.averages = averages;
    out.rawAverages = rawAverages;
    out.subspecs = subspecs;
    out.rawSubspecs = rawSubspecs;
    out.seq = seq;
    out.te = TE;
    out.tr = TR;
    out.pointsToLeftshift = pointsToLeftShift;
    out = setFov(out, 'x', rda.FoVHeight);
    out = setFov(out, 'y', rda.FoVWidth);
    out = setFov(out, 'z', rda.FoV3D);
    out = setVoxelSize(out, 'x', getFov(out, 'x')/getSizeFromDimensions(out, {'x'}));
    out = setVoxelSize(out, 'y', getFov(out, 'y')/getSizeFromDimensions(out, {'y'}));
    if(getDimension(out, 'z') ~= 0)
        out = setVoxelSize(out, 'z', getFov(out, 'z')/getSizeFromDimensions(out, {'z'}));
    else
        out = setVoxelSize(out, 'z', getFov(out, 'z'));
    end

    affine_matrix = calculateAffineMatrix(rda);
    out.affineMatrix = affine_matrix;
    out = setImageOrigin(out, rda.imageOrigin);
    out = findCoordinates(out);
    out = setFlags(out);
end

function affine_matrix = calculateAffineMatrix(rda)
    %may have to update x and y orientations
    affine_matrix_rotation = [rda.RowVector', rda.ColumnVector', rda.sliceVector', [0 0 0]'];
    affine_matrix_rotation = cat(1, affine_matrix_rotation,  [0, 0, 0,1]);
    
    affine_matrix_translate = eye(4);
    affine_matrix_translate(1,4) = rda.PositionVector(1);
    affine_matrix_translate(2,4) = rda.PositionVector(2);
    affine_matrix_translate(3,4) = rda.PositionVector(3);

    affine_matrix_scale = eye(4);
    affine_matrix_scale(1,1) = rda.PixelSpacingCol;
    affine_matrix_scale(2,2) = rda.PixelSpacingRow;
    affine_matrix_scale(3,3) = rda.PixelSpacing3D;
    affine_matrix = affine_matrix_rotation * affine_matrix_translate * affine_matrix_scale;
end

function out = setFlags(out)
    %FILLING IN THE FLAGS
    out.flags.writtentostruct = 1;
    out.flags.gotparams = 1;
    out.flags.leftshifted = 0;
    out.flags.filtered = 0;
    out.flags.zeropadded = 0;
    out.flags.freqcorrected = 0;
    out.flags.phasecorrected = 0;
    out.flags.averaged = 1;
    out.flags.addedrcvrs = 1;
    out.flags.subtracted = 0;
    out.flags.writtentotext = 0;
    out.flags.downsampled = 0;
    out.flags.spatialFT = 1;
    out.flags.spectralFT = 0;
    out.flags.isCartesian = 1;
    if out.dims.subSpecs == 0
        out.flags.isISIS = 0;
    else
        out.flags.isISIS = (out.sz(out.dims.subSpecs)==4);
    end
end

function out = findCoordinates(out)
    xCoordinate = createCoordinates(getFov(out, 'x')/2, getVoxSize(out, 'x'));
    yCoordinate = createCoordinates(getFov(out, 'y')/2, getVoxSize(out, 'y'));
    zCoordinate = createCoordinates(getFov(out, 'z')/2, getVoxSize(out, 'z'));
    out = setCoordinates(out, 'x', xCoordinate);
    out = setCoordinates(out, 'y', yCoordinate);
    out = setCoordinates(out, 'z', zCoordinate);
end
