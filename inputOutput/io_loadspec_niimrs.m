% io_loadspec_niimrs.m
% Georg Oeltzschner, Johns Hopkins University 2021
%
% USAGE:
% out = io_loadspec_niimrs(filename);
% 
% DESCRIPTION:
% Reads in MRS data stored according to the NIfTI MRS format.
% See the specification under
% https://docs.google.com/document/d/1tC4ugzGUPLoqHRGrWvOcGCuCh_Dogx_uu0cxKub0EsM/edit
% 
% io_loadspec_niimrs outputs the data in structure format, with fields 
% corresponding to time scale, fids, frequency scale, spectra, and header 
% fields containing information about the acquisition. The resulting matlab
% structure can be operated on by the other functions in this MRS toolbox.
%
% This function is currently work-in-progress and is being tested on more
% and more datasets. Please contact the FID-A developers for help with a
% particular NIfTI-MRS file that you might encounter problems with when
% using this function.
%
% Currently, this function is limited to single-voxel MRS data. It is
% planned to develop it for full compatibility with 2D and 3D multi-voxel
% and spectroscopic imaging data.
%
% DEPENDENCIES:
% This function requires the dcm2nii toolbox (Xiangrui Li) to be on the
% MATLAB path
% https://github.com/xiangruili/dicm2nii
% 
% INPUTS:
% filename   = filename of NIfTI MRS file (*.nii) to load.
%
% OUTPUTS:
% out        = Input dataset in FID-A structure format.

function out = io_loadspec_niimrs(filename)

% Read in the data using the dicm2nii toolbox
% (https://github.com/xiangruili/dicm2nii)
try
    nii = nii_tool('load', filename);
catch ME
    switch ME.identifier
        case 'MATLAB:UndefinedFunction'
            error(['Cannot find the function ''nii_tool.m''.' ...
                ' Please ensure that you have downloaded the required', ...
                ' dcm2nii toolbox (https://github.com/xiangruili/dicm2nii)', ...
                ' and added it to your MATLAB path.']);
        otherwise
            rethrow(ME);
    end
end

% Extract the header and header extensions
hdr = nii.hdr;
hdr_ext = jsondecode(nii.ext.edata_decoded);

% Extract the time-domain data
fids = nii.img;

% Extract spectrometer frequency and dwell time
f0 = hdr_ext.SpectrometerFrequency;
dt = hdr.pixdim(5);
sw = 1/dt;

% Specify dimensions
% In NIfTI MRS, the three spatial dimensions and the time dimension occupy
% fixed indices in the (maximum) 7-D array
dims.x = 1;
dims.y = 2;
dims.z = 3;
dims.t = 4;

% There are some pre-defined dimension names according to the FID-A
% convention. These dimensions may or may not be stored in the NIfTI MRS
% header, so we'll initialize them as 0.
dims.coils = 0;
dims.averages = 0;
dims.subSpecs = 0;
dims.extras = 0;

% The NIfTI MRS standard reserves the remaining 3 dimensions, which are
% then explicitly specified in the JSON header extension fields dim_5,
% dim_6 and dim_7.
if isfield(hdr_ext, 'dim_5')
    dim_number = 5;
    switch hdr_ext.dim_5
        case 'DIM_COIL'
            dims.coils      = dim_number;
        case 'DIM_DYN'
            dims.averages   = dim_number;
        case 'DIM_INDIRECT_0'
            dims.extras     = dim_number;
        case 'DIM_INDIRECT_1'
            dims.extras     = dim_number;
        case 'DIM_INDIRECT_2'
            dims.extras     = dim_number;
        case 'DIM_PHASE_CYCLE'
            dims.extras     = dim_number;
        case 'DIM_EDIT'
            dims.subSpecs   = dim_number;
        case 'DIM_MEAS'
            dims.extras     = dim_number;
        case 'DIM_USER_0'
            dims.extras     = dim_number;
        case 'DIM_USER_1'
            dims.extras     = dim_number;
        case 'DIM_USER_2'
            dims.extras     = dim_number;
        case 'DIM_ISIS'
            dims.subSpecs   = dim_number;
        otherwise
            error('Unknown dimension value specified in dim_5: %s', hdr_ext.dim_5);
    end
       
end
if isfield(hdr_ext, 'dim_6')
        dim_number = 6;
    switch hdr_ext.dim_6
        case 'DIM_COIL'
            dims.coils      = dim_number;
        case 'DIM_DYN'
            dims.averages   = dim_number;
        case 'DIM_INDIRECT_0'
            dims.extras     = dim_number;
        case 'DIM_INDIRECT_1'
            dims.extras     = dim_number;
        case 'DIM_INDIRECT_2'
            dims.extras     = dim_number;
        case 'DIM_PHASE_CYCLE'
            dims.extras     = dim_number;
        case 'DIM_EDIT'
            dims.subSpecs   = dim_number;
        case 'DIM_MEAS'
            dims.extras     = dim_number;
        case 'DIM_USER_0'
            dims.extras     = dim_number;
        case 'DIM_USER_1'
            dims.extras     = dim_number;
        case 'DIM_USER_2'
            dims.extras     = dim_number;
        case 'DIM_ISIS'
            dims.subSpecs   = dim_number;
        otherwise
            error('Unknown dimension value specified in dim_6: %s', hdr_ext.dim_6);
    end
end
if isfield(hdr_ext, 'dim_7')
        dim_number = 7;
    switch hdr_ext.dim_7
        case 'DIM_COIL'
            dims.coils      = dim_number;
        case 'DIM_DYN'
            dims.averages   = dim_number;
        case 'DIM_INDIRECT_0'
            dims.extras     = dim_number;
        case 'DIM_INDIRECT_1'
            dims.extras     = dim_number;
        case 'DIM_INDIRECT_2'
            dims.extras     = dim_number;
        case 'DIM_PHASE_CYCLE'
            dims.extras     = dim_number;
        case 'DIM_EDIT'
            dims.subSpecs   = dim_number;
        case 'DIM_MEAS'
            dims.extras     = dim_number;
        case 'DIM_USER_0'
            dims.extras     = dim_number;
        case 'DIM_USER_1'
            dims.extras     = dim_number;
        case 'DIM_USER_2'
            dims.extras     = dim_number;
        case 'DIM_ISIS'
            dims.subSpecs   = dim_number;
        otherwise
            error('Unknown dimension value specified in dim_7: %s', hdr_ext.dim_7);
    end
end

% Parse the NIfTI hdr.dim field:
allDims = hdr.dim(2:end); % all dimensions (including singletons)

% Find the number of points
nPts = allDims(dims.t);

% Find the number of averages.  'averages' will specify the current number
% of averages in the dataset as it is processed, which may be subject to
% change.  'rawAverages' will specify the original number of acquired 
% averages in the dataset, which is unchangeable.
if dims.subSpecs ~= 0
    if dims.averages ~= 0
        averages = allDims(dims.averages)*allDims(dims.subSpecs);
        rawAverages = averages;
    else
        averages = allDims(dims.subSpecs);
        rawAverages = 1;
    end
else
    if dims.averages ~= 0
        averages = allDims(dims.averages);
        rawAverages = averages;
    else
        averages = 1;
        rawAverages = 1;
    end
end

% FIND THE NUMBER OF SUBSPECS
% 'subspecs' will specify the current number of subspectra in the dataset 
% as it is processed, which may be subject to change. 'rawSubspecs' will 
% specify the original number of acquired  subspectra in the dataset, which
% is unchangeable.
if dims.subSpecs ~=0
    subspecs = allDims(dims.subSpecs);
    rawSubspecs = subspecs;
else
    subspecs = 1;
    rawSubspecs = subspecs;
end

% ORDERING THE DATA AND DIMENSIONS
% The FID-A array ordering conventions differ from the NIfTI MRS
% convention. Most importantly, FID-A is primarily tailored towards
% single-voxel data. We will start designing this function towards this
% purpose, and later adapt the formalism proposed in the csi_mod branch of
% the FID-A GitHub repository.

if allDims(1)*allDims(2)*allDims(3) == 1 % x=y=z=1
    dims.x = 0;
    dims.y = 0;
    dims.z = 0;
    fids = squeeze(fids);
    
    %Now that we've indexed the dimensions of the data array, we now need to
    %permute it so that the order of the dimensions is standardized:  we want
    %the order to be as follows:
    %   1) time domain data.
    %   2) coils.
    %   3) averages.
    %   4) subSpecs.
    %   5) extras.

    % Adjust dimension indices for the fact that we have collapsed the
    % three spatial dimensions (which we don't need for SVS data)
    sqzDims = {};
    dimsFieldNames = fieldnames(dims);
    for rr = 1:length(dimsFieldNames)
        if dims.(dimsFieldNames{rr}) ~= 0
            % Subtract 3 (x, y, z) from the dimension indices
            dims.(dimsFieldNames{rr}) = dims.(dimsFieldNames{rr}) - 3;
            sqzDims{end+1} = dimsFieldNames{rr};
        end
    end

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
        dims.t=1;dims.coils=0;dims.averages=0;dims.subSpecs=0;dims.extras=0;
    end
    
    %Now get the size of the data array:
    sz=size(fids);
    
    %Compared to NIfTI MRS, FID-A needs the conjugate
    fids = conj(fids);
    
    %Now take fft of time domain to get fid:
    specs=fftshift(ifft(fids,[],dims.t),dims.t);

end



% Fill in additional FID-A format variables
% Nucleus (new field)
out.nucleus = hdr_ext.ResonantNucleus;
% Calculate B0 from spectrometer frequency depending on nucleus
% Gamma from Wikipedia article "Gyromagnetic ratio" (3 signif. digits)
for rr = 1:length(out.nucleus)
    switch out.nucleus{rr}
        case '1H'
            gamma = 42.577;
        case '2H'
            gamma = 6.536;
        case '3HE'
            gamma = -32.434;
        case '7LI'
            gamma = 16.546;
        case '13C'
            gamma = 10.708;
        case '19F'
            gamma = 40.052;
        case '23NA'
            gamma = 11.262;
        case '31P'
            gamma = 17.235;
        case '129XE'
            gamma = -11.777;
    end
    Bo(rr) = f0(rr) ./ gamma;
end

% Calculate t and ppm arrays using the calculated parameters:
f   =[(-sw/2) + (sw/(2*nPts)) : sw/(nPts) : (sw/2) - (sw/(2*nPts))];
ppm = -f / (Bo(1)*42.577);
ppm = ppm + 4.65;
t   = [0 : dt : (nPts-1)*dt];


% MANDATORY FIELDS
% Data & dimensions
out.fids = fids;
out.specs = specs;
out.sz = sz;
out.dims = dims;
out.Bo = Bo;
out.averages = averages;
out.rawAverages = rawAverages;
out.subspecs = subspecs;
out.rawSubspecs = rawSubspecs;

% Echo/repetition time
out.te = hdr_ext.EchoTime;
out.tr = hdr_ext.RepetitionTime;

% time and frequency axis
out.t   = t;
out.ppm = ppm;

% Dwell time & spectral width & field strength
out.spectralwidth = sw;
out.dwelltime = dt;
out.txfrq  = f0 * 10^6;
out.date = '';

% NIfTI-MRS-SPECIFIC FIELDS
% Save the NIfTI header
out.nii_mrs.hdr = hdr;
% Save the header extension
out.nii_mrs.hdr_ext = hdr_ext;
if isfield(hdr_ext, 'SequenceName')
    out.seq = hdr_ext.SequenceName;
end


%FILLING IN THE FLAGS
out.flags.writtentostruct=1;
out.flags.gotparams=1;
out.flags.leftshifted=0;
out.flags.filtered=0;
out.flags.zeropadded=0;
out.flags.freqcorrected=0;
out.flags.phasecorrected=0;
if out.dims.averages==0
    out.flags.averaged=1;
else
    out.flags.averaged=0;
end
if out.dims.coils==0
    out.flags.addedrcvrs=1;
else
    out.flags.addedrcvrs=0;
end
out.flags.subtracted=0;
out.flags.writtentotext=0;
out.flags.downsampled=0;
if out.dims.subSpecs==0
    out.flags.isFourSteps=0;
else
    out.flags.isFourSteps=(out.sz(out.dims.subSpecs)==4);
end

end

