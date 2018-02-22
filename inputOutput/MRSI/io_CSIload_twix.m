%io_CSIload_twix.m
%Brenden Kadota, Jamie Near, Sunnybrook 2021.
%
% USAGE:
% out=io_CSIload_twix('filename');
% out=io_CSIload_twix(scan_number);
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

function out=io_CSIload_twix(filename)
arguments
    filename (1,:) {mustBeFile}
end


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
    twix_obj=mapVBVD(char(filename));
end

if isstruct(twix_obj)
    disp('single RAID file detected.');
elseif iscell(twix_obj)
    disp('multi RAID file detected.');
    RaidLength=length(twix_obj);
    %assume that the data of interest is in the last element of the cell.
    twix_obj=twix_obj{RaidLength};
end

%get CSI fids
dOut.data=twix_obj.image();
%Squeeze the data to remove singleton dims
fids=squeeze(dOut.data);

%get version number and dims
version=twix_obj.image.softwareVersion;
sqzDims=twix_obj.image.sqzDims;

%find out what sequence, the data were acquired with.  If this is a
%multi-raid file, then the header may contain multiple instances of
%'tSequenceFileName' for different scans (including a pre-scan).
%Therefore, if multi-raid file, we will need to do a bit of extra digging 
%to find the correct sequence name.  
sequence=twix_obj.hdr.Config.SequenceFileName;  

%Try to find out what sequnece this is:
isSiemens=(contains(sequence,'csi_se') ||... %Or the Siemens CSI PRESS sequence?
            contains(sequence,'csi_st'));




%Make a pulse sequence identifier for the header (out.seq);
seq=sequence;

%Find the magnetic field strength:
Bo=twix_obj.hdr.Dicom.flMagneticFieldStrength;


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
dims.x = find(strcmp(sqzDims,'Seg') | strcmp(sqzDims, 'Phs'));
if ~isempty(dims.x)
    %remove the coils dimension from the dimsToIndex vector
    dimsToIndex=dimsToIndex(dimsToIndex~=dims.x);
else
    dims.x=0;
end

%k space coordinates in the y direction
dims.y = find(strcmp(sqzDims,'Lin') );
if ~isempty(dims.y)
    %remove the coils dimension from the dimsToIndex vector
    dimsToIndex=dimsToIndex(dimsToIndex~=dims.y);
else
    dims.y=0;
end

%k space coordinates in the z direction 
dims.z = find(strcmp(sqzDims,'Sli'));
if ~isempty(dims.z)
    %remove the coils dimension from the dimsToIndex vector
    dimsToIndex=dimsToIndex(dimsToIndex~=dims.z);
else
    dims.z=0;
end

if ~isempty(dimsToIndex)
    
    dims.subSpecs=dimsToIndex(1);
    
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
%   5) csi_x.
%   6) csi_y.
%   7) extras.
dimsArray = [dims.t dims.coils, dims.averages dims.subSpecs dims.x dims.y dims.extras];
leftover = [];
counter = 1;
fields = ["t", "coils", "averages", "subSpecs", "x", "y", "extras"];
dimsCell = strings(1,1);
for i = 1:numel(dimsArray)
    if(dimsArray(i) ~= 0)
        leftover(counter) = dimsArray(i);
        dimsCell(counter) = fields(i);
        counter = counter + 1;
    end
end
fids = permute(fids, leftover);
for i = 1:size(leftover, 2)
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
fovZ = twix_obj.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness;

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

%Get TxFrq
%txfrq=twix_obj.hdr.Meas.Frequency;
if(isSiemens)
    leftshift = twix_obj.image.freeParam(1);
else
    leftshift = twix_obj.image.freeParam(1);
end

txfrq=twix_obj.hdr.Meas.Frequency;

%date = getfield(regexp(twix_obj.hdr.MeasYaps.tReferenceImage0, ...
%'^".*\.(?<DATE>\d{8})\d*"$', 'names'), 'DATE');  %Franck Lamberton

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
out.fovX = fovX;    %[mm]
out.fovY = fovY;    %[mm]
out.fovZ = fovZ;    %[mm]
out.deltaX = fovX/size(fids, dims.x);
out.deltaY = fovY/size(fids, dims.y);
if dims.z == 0
    out.deltaZ = fovZ;
else
    out.deltaZ = fovZ/size(fids, dims.z);
end
out.averages = averages;
out.rawAverages = rawAverages;
out.subspecs = subspecs;
out.rawSubspecs = rawSubspecs;

out.pointsToLeftshift=leftshift;

out.imageOrigin = zeros(1,3);
if(isfield(twix_obj.hdr.Config, 'VoI_Position_Sag'))
    fields = ["VoI_Position_Sag", "VoI_Position_Cor", "VoI_Position_Tra"];
elseif(isfield(twix_obj.hdr.Config, 'Voi_Position_Sag'))
    fields = ["Voi_Position_Sag", "Voi_Position_Cor", "Voi_Position_Tra"];
else
    fields = [];
end

for i = 1:length(fields)
    if(~isempty(twix_obj.hdr.Config.(fields(i))))
        out.imageOrigin(i) = twix_obj.hdr.Config.(fields(i));
    end
end

%initalize to zero if empty
z_vect(1) = -initalize_zero_if_emtpy(twix_obj.hdr.Meas.VoI_Normal_Sag);
z_vect(2) = -initalize_zero_if_emtpy(twix_obj.hdr.Meas.VoI_Normal_Cor);
z_vect(3) = initalize_one_if_emtpy(twix_obj.hdr.Meas.VoI_Normal_Tra);

theta = initalize_zero_if_emtpy(twix_obj.hdr.Meas.VoI_InPlaneRotAngle);
%get rotation matrix from vector and rotation angle
rot_matrix = rot_vector_to_rot_matrix(z_vect(1), z_vect(2), z_vect(3), theta);
%create affine matrix from rotation matrix
rot_matrix(4,4) = 1;

%find y vector perpendicular to z with x equal to zero
y_vect = [0, 1, -z_vect(2)/z_vect(3)];
%rotate y vector around z using rotation matrix
y_vect = y_vect/norm(y_vect);
y_vect = rot_matrix*[y_vect'; 1];
y_vect = y_vect(1:3);
%get x vector from taking the cross product of y and z. Ensures x y z are
%perpendicular
x_vect = cross(y_vect, z_vect);

%get affine rotation matrix
affine_rotation_matrix = [x_vect', y_vect, z_vect', [0,0,0]'; 0 0 0 1];

%get scaling matrix
affine_scale_matrix = eye(4);
affine_scale_matrix(1,1) = out.deltaX;
affine_scale_matrix(2,2) = -out.deltaY;
affine_scale_matrix(3,3) = out.deltaZ;

%get translate matrix
affine_translate_matrix = eye(4);
affine_translate_matrix(1,4) = -out.imageOrigin(1) - out.fovX/2;
affine_translate_matrix(2,4) = -out.imageOrigin(2) + out.fovY/2;
affine_translate_matrix(3,4) = out.imageOrigin(3) - out.fovZ/2;

%ORDER MATTERS HERE!!! First we scale image coordinates to be voxel sizes.
%Then we shift the voxels so the first voxel starts at smallest x,y,z
%coordinate. Then we rotate the voxels using the rotation matrix
out.affine_matrix = affine_rotation_matrix*affine_translate_matrix*affine_scale_matrix;

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
out.flags.isFourSteps = 0;
end


%convert rotation around vector to rotation matrix. Formula found from
%wikipedia:
%https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle.
function rotation_matrix = rot_vector_to_rot_matrix(x, y, z, theta)
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

function out = initalize_zero_if_emtpy(value)
    if(isempty(value))
        out = 0;
    else
        out = value;
    end
end

function out = initalize_one_if_emtpy(value)
    if(isempty(value))
        out = 1;
    else
        out = value;
    end
end