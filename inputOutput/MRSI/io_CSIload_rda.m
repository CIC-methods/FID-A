%io_CSIload_rda.m
%Brenden Kadota, SunnyBrook Hosptial 2021.
%
% USAGE:
% [out]=io_CSIload_rda(rda_filename);
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
            eval(['rda.' , variable , ' = value; ']);
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
            eval(['rda.' , variable , ' = str2num(value); ']);
        case { 'PatientWeight' , 'TR' , 'TE' , 'TM' , 'DwellTime' , 'NumberOfAverages' , 'MRFrequency' , 'MagneticFieldStrength' , 'FlipAngle' , ...
                     'SliceThickness' ,  'FoVHeight' , 'FoVWidth', 'FoV3D' , 'PercentOfRectFoV' , 'PixelSpacingRow' , 'PixelSpacingCol', 'PixelSpacing3D'}
            %Floats 
            eval(['rda.' , variable , ' = str2num(value); ']);
        case {'SoftwareVersion[0]' }
            rda.software_version = value;
        case {'CSIMatrixSize[0]' }
            rda.CSIMatrix_Size(1) = str2num(value);    
        case {'CSIMatrixSize[1]' }
            rda.CSIMatrix_Size(2) = str2num(value);    
        case {'CSIMatrixSize[2]' }
            rda.CSIMatrix_Size(3) = str2num(value);   
        case {'CSIMatrixSizeOfScan[0]'}
            rda.CSIMatrixSizeOfScan(1) = str2num(value);
        case {'CSIMatrixSizeOfScan[1]'}
            rda.CSIMatrixSizeOfScan(2) = str2num(value);
        case {'CSIMatrixSizeOfScan[2]'}
            rda.CSIMatrixSizeOfScan(3) = str2num(value);
        case {'PositionVector[0]' }
            rda.PositionVector(1) = str2num(value);    
        case {'PositionVector[1]' }
            rda.PositionVector(2) = str2num(value);     
        case {'PositionVector[2]' }
            rda.PositionVector(3) = str2num(value);    
        case {'RowVector[0]' }
            rda.RowVector(1) = str2num(value);    
        case {'RowVector[1]' }
            rda.RowVector(2) = str2num(value);
        case {'RowVector[2]' }
            rda.RowVector(3) = str2num(value);
        case {'ColumnVector[0]' }
            rda.ColumnVector(1) = str2num(value);
        case {'ColumnVector[1]' }
            rda.ColumnVector(2) = str2num(value);
        case {'ColumnVector[2]' }
            rda.ColumnVector(3) = str2num(value);
        case {'VOIPositionSag' }
            rda.imageOrigin(1) = str2num(value);
        case {'VOIPositionCor' }
            rda.imageOrigin(2) = str2num(value);
        case {'VOIPositionTra' }
            rda.imageOrigin(3) = str2num(value);
        case {'VOINormalSag' }
            rda.sliceVector(1) = str2num(value);
        case {'VOINormalCor' }
            rda.sliceVector(2) = str2num(value);
        case {'VOINormalTra' }
            rda.sliceVector(3) = str2num(value);
                
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
 fids = flip(fids, 2);
 fids = conj(fids);

% get the spectrum from the fid

% make calculations for the output mrs structure
sz = size(fids);
dwelltime = rda.DwellTime/1000000;
spectralwidth=1/dwelltime;
txfrq = rda.MRFrequency*1000000;
dims.t = 1;
dims.subSpecs = 0;
dims.x = 2;
dims.y = 3;
dims.z = 0;
dims.coils = 0;
dims.averages = 0;
Bo = rda.MagneticFieldStrength;
rawAverages = rda.NumberOfAverages;
averages = 1;
subspecs =1;
rawSubspecs = 'na';
date = rda.StudyDate;
seq = rda.SequenceDescription;
TE = rda.TE;
TR = rda.TR;
pointsToLeftShift = 'N/A';

%Calculate t and ppm arrays using the calculated parameters:
f=[(-spectralwidth/2)+(spectralwidth/(2*sz(1))):spectralwidth/(sz(1)):(spectralwidth/2)-(spectralwidth/(2*sz(1)))];
ppm=-f/(Bo*42.577);
ppm=ppm+4.6082;

t=[0:dwelltime:(sz(1)-1)*dwelltime];

%FILLING IN DATA STRUCTURE
out.fids=fids;
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
out.te=TE;
out.tr=TR;
out.pointsToLeftshift=pointsToLeftShift;
out.fovY = rda.FoVWidth;    
out.fovX = rda.FoVHeight;   
out.fovZ = rda.FoV3D;
out.deltaX = out.fovX/size(fids, dims.x);
out.deltaY = out.fovY/size(fids, dims.y);
out.deltaZ = out.fovZ;



%may have to update x and y orientations
affine_matrix_rotation = [rda.RowVector', rda.ColumnVector', rda.sliceVector', [0 0 0]'];
affine_matrix_rotation = cat(1, affine_matrix_rotation,  [0, 0, 0,1]);

affine_matrix_translate = eye(4);
affine_matrix_translate(1,4) = rda.PositionVector(1);
affine_matrix_translate(2,4) = -rda.PositionVector(2);
affine_matrix_translate(3,4) = rda.PositionVector(3);

affine_matrix_scale = eye(4);
affine_matrix_scale(1,1) = rda.PixelSpacingCol; 
affine_matrix_scale(2,2) = -rda.PixelSpacingRow;
affine_matrix_scale(3,3) = rda.PixelSpacing3D;
out.affine_matrix = affine_matrix_rotation * affine_matrix_translate * affine_matrix_scale;

out.imageOrigin = rda.imageOrigin;
out.xCoordinates = (-out.fovX/2 + out.deltaX/2:out.deltaX:out.fovX/2 - out.deltaX/2) + out.imageOrigin(1);    
out.yCoordinates = (-out.fovY/2 + out.deltaY/2:out.deltaY:out.fovY/2 - out.deltaY/2) + out.imageOrigin(2);   
out.zCoordinates = (-out.fovZ/2 + out.deltaZ/2:out.deltaZ:out.fovZ/2 - out.deltaZ/2) + out.imageOrigin(3);


%FILLING IN THE FLAGS
out.flags.writtentostruct=1;
out.flags.gotparams=1;
out.flags.leftshifted=0;
out.flags.filtered=0;
out.flags.zeropadded=0;
out.flags.freqcorrected=0;
out.flags.phasecorrected=0;
out.flags.averaged=1;
out.flags.addedrcvrs=1;
out.flags.subtracted=0;
out.flags.writtentotext=0;
out.flags.downsampled=0;
out.flags.spatialFT = 1;
out.flags.spectralFT = 0;
if out.dims.subSpecs==0
    out.flags.isISIS=0;
else
    out.flags.isISIS=(out.sz(out.dims.subSpecs)==4);
end
      
end
