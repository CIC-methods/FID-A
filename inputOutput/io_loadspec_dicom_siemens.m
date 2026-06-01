%io_loadspec_dicom_siemens.m
%Jamie Near, Sunnybrook Research Institute, 2026.
%
% USAGE:
% [out,info]=io_loadspec_dicom_siemens(dicomFolder);
%
% DESCRIPTION:
% Reads in an MRS dataset in Siemens dicom format.  The input is a folder 
% of Siemens dicom files, and the output is a single FID-A structure.  The
% resulting matlab structure can be operated on by the other functions in
% the FID-A toolkit. 

% Supports:
%   - Siemens XA/XA20/XA30/XA50/XA60 MR Spectroscopy Storage DICOMs
%   - Siemens VE11/E11 private CSA spectroscopy DICOMs
%
% INPUTS:
% dicomFolder   = Folder of Siemens dicom files to load
%
% OUTPUTS: 
% out           = Input dataset in FID-A structure format
% info          = Complete dicom header information
%

function [raw,info0] = io_loadspec_dicom_siemens(folder)

%% ------------------------------------------------------------
%  Get files
%% ------------------------------------------------------------

files = dir(folder);
files = files(~[files.isdir]);

if isempty(files)
    error('No files found in folder.');
end

% Remove hidden files
keep = true(size(files));
for f = 1:length(files)
    if startsWith(files(f).name,'.')
        keep(f) = false;
    end
end
files = files(keep);
nFil=length(files);

%% ------------------------------------------------------------
%  Read headers and sort by InstanceNumber
%% ------------------------------------------------------------

instanceNumbers = zeros(length(files),1);
infos = cell(length(files),1);

fprintf('Reading DICOM headers...\n');

for f = 1:nFil;

    fname = fullfile(folder,files(f).name);

    try
        info = dicominfo(fname);
    catch
        warning('Skipping non-DICOM file: %s',files(f).name);
        continue;
    end

    infos{f} = info;

    if isfield(info,'InstanceNumber')
        instanceNumbers(f) = double(info.InstanceNumber);
    else
        instanceNumbers(f) = f;
    end
end

% Remove empties
valid = ~cellfun(@isempty,infos);
infos = infos(valid);
files = files(valid);
instanceNumbers = instanceNumbers(valid);

% Sort by InstanceNumber
[~,ord] = sort(instanceNumbers);
files = files(ord);
infos = infos(ord);

nAvg = length(files);

fprintf('Found %d spectroscopy files.\n',nAvg);

%% ------------------------------------------------------------
%  Determine format
%% ------------------------------------------------------------

info0 = infos{1};

isXA = isfield(info0,'SpectroscopyData');
isVE = isfield(info0,'Private_7fe1_1010');

if ~isXA && ~isVE
    error('Unsupported Siemens DICOM spectroscopy format.');
end

%XA FORMAT FIRST
if isXA

    fprintf('Detected XA-style spectroscopy DICOM.\n');

    % Get the CSA information directly from the dicom file (not from dicomInfo)
    fname = fullfile(folder,files(1).name);
    csa = extractPhoenixFromDicomFile(fname);

    % ---------- read data ----------
    for a = 1:nAvg
        info = infos{a};
        data = double(info.SpectroscopyData);
        data = data(1:2:end) + 1i*data(2:2:end);
        data = data(:);
        fids(:,a) = data;
    end


% NOW VE11 / E11 FORMAT
elseif isVE

    fprintf('Detected VE11/E11 private CSA spectroscopy DICOM.\n');

    % Read CSA header 
    csa = extractPhoenixfromDicomInfo(info0);

    % Read FIDs
    for a = 1:nAvg
        info = infos{a};
        rawbytes = info.Private_7fe1_1010;
        tmp = typecast(rawbytes,'single');
        data = double(tmp(1:2:end)) + 1i*double(tmp(2:2:end));
        data = data(:);
        fids(:,a) = data;
    end
end

%Make sz variable:
sz=size(fids);

%Get vector size:
nPts = getPhoenixValue(csa,'lVectorSize',true);

%Get dwell time:
dwelltime = getPhoenixValue(csa,'alDwellTime[0]',true);
dwelltime = dwelltime*1e-9;
spectralwidth = 1/dwelltime;

%Get frequency
txfrq = getPhoenixValue(csa,'lFrequency',true);

%Get B0
Bo = getPhoenixValue(csa,'flNominalB0',true);

%Get TE and TR
te = getPhoenixValue(csa,'alTE[0]',true);
te = te * 1e-3;
tr = getPhoenixValue(csa,'alTR[0]',true);
tr = tr * 1e-3;

% Get date
date = infos{1}.PerformedProcedureStepStartDate;

%Find out what sequence this is:
sequence = getPhoenixValue(csa,'tSequenceFileName',false);

%Check if this is edited data:
if contains(sequence,'eja_svs_mpress') || contains(sequence,'special')
    fids=reshape(fids,[sz(1),sz(2)/2,2]);
    averages = nFil/2;
    rawAverages = nFil/2;
    subspecs = 2;
    rawSubspecs = 2;

    %Make the dims structure:
    dims.t = 1;
    dims.coils=0;
    dims.averages=2;
    dims.subSpecs=3;
    dims.extras=0;

else
    averages = nFil;
    rawAverages= nFil;
    subspecs = 1;
    rawSubspecs = 1;

    %Make the dims structure:
    dims.t = 1;
    dims.coils=0;
    dims.averages=2;
    dims.subSpecs=0;
    dims.extras=0;
end

%% ------------------------------------------------------------
%  Spectral domain
%% ------------------------------------------------------------

%Calculate ppm and specs array using the calculated parameters:
f=[(-spectralwidth/2)+(spectralwidth/(2*sz(1))):spectralwidth/(sz(1)):(spectralwidth/2)-(spectralwidth/(2*sz(1)))];
ppm=-f/(Bo*42.577);
ppm=ppm+4.6082;

specs = fftshift(ifft(fids,[],1),1);


%% ------------------------------------------------------------
%  Build FID-A structure
%% ------------------------------------------------------------

raw = struct();

raw.fids = fids;
raw.specs = specs;
raw.sz = size(fids);
raw.ppm = ppm;
raw.t = (0:nPts-1) * dwelltime;
raw.spectralwidth = spectralwidth;
raw.dwelltime = dwelltime;
raw.txfrq = txfrq;
raw.date=date;
raw.dims=dims;
raw.Bo = Bo;
raw.averages = nAvg;
raw.rawAverages = nAvg;
raw.subspecs = 1;
raw.rawSubspecs = 1;
raw.seq = sequence;
raw.te = te;
raw.tr = tr;
raw.pointsToLeftshift = 0;
raw.nucleus = '1H';
raw.flags = struct();
raw.flags.writtentostruct = 1;
raw.flags.gotparams = 1;
raw.flags.leftshifted = 0;
raw.flags.filtered = 0;
raw.flags.zeropadded = 0;
raw.flags.freqcorrected = 0;
raw.flags.phasecorrected = 0;
raw.flags.averaged = 0;
raw.flags.addedrcvrs = 1;
raw.flags.subtracted = 0;
raw.flags.writtentotext = 0;
raw.flags.downsampled = 0;
raw.flags.isISIS = 0;

end

%% =====================================================================
%  Helper functions
%% =====================================================================

function txt = extractPhoenixfromDicomInfo(info)

fields = fieldnames(info);
txt = '';
for f = 1:length(fields)
    fld = fields{f};
    try
        val = info.(fld);

        % Only inspect potentially useful private fields
        if contains(lower(fld),'private') || ...
                contains(lower(fld),'csa')

            % Numeric arrays -> char
            if isnumeric(val) || islogical(val)
                try
                    val = char(val');
                catch
                    continue;
                end
            end

            % Cell arrays
            if iscell(val)
                continue;
            end

            % Convert to char safely
            val = char(val(:)');

            % Keep only fields containing protocol text
            if contains(val,'ASCCONV') || ...
                    contains(val,'alTE') || ...
                    contains(val,'tSequenceFileName')

                txt = [txt val]; %#ok<AGROW>

            end
        end

    catch
    end
end

end

%% ---------------------------------------------------------------------

function val = parse_csa_numeric(txt,pattern)

val = [];

expr = [pattern '.*?([-+]?[0-9]*\.?[0-9]+)'];

tok = regexp(txt,expr,'tokens','once');

if isempty(tok)
    return;
end

val = str2double(tok{1});

end

%% ---------------------------------------------------------------------

function txt = extractPhoenixFromDicomFile(fname)

fid = fopen(fname,'r');

raw = fread(fid,'uint8=>char')';

fclose(fid);

i1 = strfind(raw,'ASCCONV BEGIN');
i2 = strfind(raw,'ASCCONV END');

if isempty(i1) || isempty(i2)
    txt = '';
    return
end

txt = raw(i1(1):i2(1)+length('ASCCONV END'));

end 

%% --------------------------------------------------------------------
function value = getPhoenixValue(txt,param,isNumeric)

value = [];

% Find parameter
idx = strfind(txt,param);

if isempty(idx)
    return
end

% Find '=' after parameter
eqIdx = strfind(txt(idx:end),'=');

if isempty(eqIdx)
    return
end

eqIdx = idx + eqIdx(1) - 1;

% Find next newline
nlIdx = regexp(txt(eqIdx:end),'[\r\n]','once');

if isempty(nlIdx)
    return
end

nlIdx = eqIdx + nlIdx - 2;

% Extract text between '=' and newline
value = strtrim(txt(eqIdx+1:nlIdx));

if isNumeric
    % Convert to numeric
    value = str2double(value);
end

end

%% -----------------------------------------------------------------------