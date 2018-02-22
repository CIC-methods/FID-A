function twix_obj = mapVBVD(filename,varargin)

%  Reads Siemens raw .dat file from VB/VD MRI raw data.
%  requires twix_map_obj.m & read_twix_hdr.m
%
%
%  Author: Philipp Ehses (philipp.ehses@tuebingen.mpg.de)
%  
% 
%  Philipp Ehses 11.02.07, original version
%  [..]
%  Philipp Ehses 22.03.11, port to VD11
%  Felix Breuer  31.03.11, added support for noise & ref scans, speed fixes
%  Philipp Ehses 19.08.11, complete reorganization of code, added
%                          siemens_data_obj class to improve readability
%  Wolf Blecher  15.05.12, readout of slice position parameters for VB Data sets
%  Wolf Blecher  11.06.12, added distinction between PATREF and PATREF PHASCOR
%  Philipp Ehses 02.08.12, again massive code reorganization. The new class
%                          twix_map_obj.m now stores the memory position of
%                          each dataset (coils are included - size: NCol*NCha)
%                          The actual data is not read until it is demanded
%                          by a "data_obj()" call. This makes it possible
%                          to selectively read in only parts of the data
%                          (e.g. to preserve memory).
%  Philipp Ehses 27.08.12, speed optimizations (avoiding of .-subsref calls)
%  07.09.12 Thanks to Stephen Yutzy for implementing support for raw data
%           correction (currently only supported for VB software version).
%  15.01.13 Thanks to Zhitao Li for proper handling of SYNCDATA.
%  Philipp Ehses 28.08.13, added support for VD13 multi-raid files
%  Michael Völker May-Aug 15 * Better error tolerance with incomplete files.
%                          * Swapped out parsing loop into a separate function
%                            without access to twix object (no thousands of subsref calls).
%                          * For parsing, use an as-minimalistic-as-possible loop
%                            to gather all mdhs in binary form. They are all stored
%                            in one array "mdh_blob".
%                          * Translation of mdhs from binary to struct is vectorized 
%                            and almost instant. It's done in evalMDH(), replacing
%                            evalMDHvb() and evalMDHvd(), no more freads inside!
%                          * vectorized readMDH(), quasi-instant
%                          * When parsing, actually read the entire file without jumps.
%                            This is weirdly faster than fseek(), plus the entire file is
%                            kept in the file system cache, if possible. Next read
%                            is therefore faster, too.
%                          * => Parsing speed improved by factor 3...7 or so
%                          * Speed increase for reading data, esp. when slicing,
%                            os-removal or reflected lines. Also for random acquisitions.
% Jonas Bause,    18.11.16   receiver phase for ramp-sampling fixed, now takes into account  
%   Chris Mirkes & PE        offcenter shifts in readout direction
% 
% 
% Input:
% 
% filename or simply meas. id, e.g. mapVBVD(122) (if file is in same path)
% optional arguments (see below)
% 
% Output: twix_obj structure with elements (if available):
%     .image:         image scan
%     .noise:         for noise scan
%     .phasecor:      phase correction scan
%     .phasestab:     phase stabilization scan
%     .phasestabRef0: phase stab. ref. (MDH_REFPHASESTABSCAN && !MDH_PHASESTABSCAN)
%     .phasestabRef1: phase stab. ref. (MDH_REFPHASESTABSCAN &&  MDH_PHASESTABSCAN)
%     .refscan:       parallel imaging reference scan
%     .refscanPC:     phase correction scan for reference data
%     .refscanPS:     phase stabilization scan for reference data
%     .refscanPSRef0: phase stab. ref scan for reference data
%     .refscanPSRef1: phase stab. ref scan for reference data
%     .RTfeedback:    realtime feedback data
%     .vop:           vop rf data
% 
% 
% The raw data can be obtained by calling e.g. twix_obj.image() or for
% squeezed data twix_obj.image{''} (the '' are needed due to a limitation
% of matlab's overloading capabilities).
% Slicing is supported as well, e.g. twix_obj.image(:,:,1,:) will return
% only the first line of the full data set (all later dimensions are
% squeezed into one). Thus, slicing of the "memory-mapped" data objects
% works exactly the same as regular matlab array slicing - with one
% exception:
% The keyword 'end' is not supported.
% Overloading of the '()' and '{}' operators works by overloading matlab's
% built-in 'subsref' function. Matlab calls subsref whenever the operators
% '()', '{}', or '.' are called. In the latter case, the overloaded subsref
% just calls the built-in subsref function since we don't want to change
% the behaviour for '.'-calls. However, this has one minor consequence:
% There's no way (afaik) to know whether the original call was terminated
% with a semicolon. Thus, a call to e.g. twix_obj.image.NLin will produce
% no output with or without semicolon termination. 'a = twix_obj.image.NLin'
% will however produce the expected result.
% 
% 
% Order of raw data:
%  1) Columns
%  2) Channels/Coils
%  3) Lines
%  4) Partitions
%  5) Slices
%  6) Averages
%  7) (Cardiac-) Phases
%  8) Contrasts/Echoes
%  9) Measurements
% 10) Sets
% 11) Segments
% 12) Ida
% 13) Idb
% 14) Idc
% 15) Idd
% 16) Ide
% 
% 
% Optional parameters/flags:
% 
% removeOS:          removes oversampling (factor 2) in read direction
% doAverage:         performs average (resulting avg-dim has thus size 1)
% ignoreSeg:         ignores segment mdh index (works basically the same as
%                    the average flag)
% rampSampRegrid     optional on-the-fly regridding of ramp-sampled readout
% doRawDataCorrect:  enables raw data correction if used in the acquisition
%                    (only works for VB atm)
% 
% These flags can also be set/unset later, e.g "twix_obj.image.flagRemoveOS = 1"
% 
% 
% Examples:
%   twix_obj = mapVBVD(measID);
% 
%   % return all image-data:
%   image_data = twix_obj.image();
%   % return all image-data with all singular dimensions removed/squeezed:
%   image_data = twix_obj.image{''}; % '' necessary due to a matlab limitation
%   % return only data for line numbers 1 and 5; all dims higher than 4 are
%   % grouped into dim 5):
%   image_data = twix_obj.image(:,:,[1 5],:,:);
%   % return only data for coil channels 2 to 6; all dims higher than 4 are
%   % grouped into dim 5); but work with the squeezed data order
%   % => use '{}' instead of '()':
%   image_data = twix_obj.image{:,2:6,:,:,:};
% 
%   So basically it works like regular matlab array slicing (but 'end' is
%   not supported; note that there are still a few bugs with array slicing).
% 
%   % NEW: unsorted raw data (in acq. order):
%   image_data = twix_obj.image.unsorted(); % no slicing supported atm
% 
% 
% Suppress silly editor warnings in the entire file, barking about
% unused variables:
%#ok<*NASGU>

if ~exist('filename','var') || isempty(filename)
    info = 'Please select binary file to read';
    [fname,pathname]=uigetfile('*.dat',info);
    if isempty(pathname)
        return
    end
    filename=[pathname fname];
else
    if ischar(filename)
        % assume that complete path is given
        if  ~strcmpi(filename(end-3:end),'.dat');
            filename=[filename '.dat'];   %% adds filetype ending to file
        end
    else
        % filename not a string, so assume that it is the MeasID
        measID   = filename;
        filelist = dir('*.dat');
        filesfound = 0;
        for k=1:numel(filelist)
            if regexp(filelist(k).name,['^meas_MID0*' num2str(measID) '_.*\.dat'])==1
                if filesfound == 0
                    filename = filelist(k).name;
                end
                filesfound = filesfound+1;
            end
        end
        if filesfound == 0
            error(['File with meas. id ' num2str(measID) ' not found.']);
        elseif filesfound > 1
            disp(['Multiple files with meas. id ' num2str(measID) ' found. Choosing first occurence.']);
        end
    end
end

% add absolute path, when no path is given
[pathstr, name, ext] = fileparts(filename);

if isempty(pathstr)
    pathstr  = pwd;
    filename = fullfile(pathstr, [name ext]);
end

%%%%% Parse varargin %%%%%
arg.bReadImaScan    = true;
arg.bReadNoiseScan  = true;
arg.bReadPCScan     = true;
arg.bReadRefScan    = true;
arg.bReadRefPCScan  = true;
arg.bReadRTfeedback = true;
arg.bReadPhaseStab  = true;
arg.bReadHeader     = true;

k=1;
while k <= numel(varargin)

    if ~ischar(varargin{k})
        error('string expected');
    end

    switch lower(varargin{k})
        case {'readheader','readhdr','header','hdr'}
            if numel(varargin) > k && ~ischar(varargin{k+1})
                arg.bReadHeader = logical(varargin{k+1});
                k = k+2;
            else
                arg.bReadHeader = true;
                k = k+1;
            end
        case {'removeos','rmos'}
            if numel(varargin) > k && ~ischar(varargin{k+1})
                arg.removeOS = logical(varargin{k+1});
                k = k+2;
            else
                arg.removeOS = true;
                k = k+1;
            end
        case {'doaverage','doave','ave','average'}
            if numel(varargin) > k && ~ischar(varargin{k+1})
                arg.doAverage = logical(varargin{k+1});
                k = k+2;
            else
                arg.doAverage = true;
                k = k+1;
            end
        case {'averagereps','averagerepetitions'}
            if numel(varargin) > k && ~ischar(varargin{k+1})
                arg.averageReps = logical(varargin{k+1});
                k = k+2;
            else
                arg.averageReps = true;
                k = k+1;
            end
        case {'averagesets'}
            if numel(varargin) > k && ~ischar(varargin{k+1})
                arg.averageSets = logical(varargin{k+1});
                k = k+2;
            else
                arg.averageSets = true;
                k = k+1;
            end
        case {'ignseg','ignsegments','ignoreseg','ignoresegments'}
            if numel(varargin) > k && ~ischar(varargin{k+1})
                arg.ignoreSeg = logical(varargin{k+1});
                k = k+2;
            else
                arg.ignoreSeg = true;
                k = k+1;
            end
        case {'rampsampregrid','regrid'}
            if numel(varargin) > k && ~ischar(varargin{k+1})
                arg.rampSampRegrid = logical(varargin{k+1});
                k = k+2;
            else
                arg.rampSampRegrid = true;
                k = k+1;
            end
        case {'rawdatacorrect','dorawdatacorrect'}
            if numel(varargin) > k && ~ischar(varargin{k+1})
                arg.doRawDataCorrect = logical(varargin{k+1});
                k = k+2;
            else
                arg.doRawDataCorrect = true;
                k = k+1;
            end
        otherwise
            error('Argument not recognized.');
    end
end
clear varargin
%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
fid = fopen(filename, 'r', 'ieee-le');

% get file size
fseek(fid,0,'eof');
fileSize = ftell(fid);

% start of actual measurement data (sans header)
fseek(fid,0,'bof');

firstInt  = fread(fid,1,'uint32');
secondInt = fread(fid,1,'uint32');

% lazy software version check (VB or VD?)
if and(firstInt < 10000, secondInt <= 64)
    version = 'vd';
    disp('Software version: VD (!?)');

    % number of different scans in file stored in 2nd in
    NScans = secondInt;
    measID = fread(fid,1,'uint32');
    fileID = fread(fid,1,'uint32');
    measOffset = cell(1, NScans);
    measLength = cell(1, NScans);
    for k=1:NScans
        measOffset{k} = fread(fid,1,'uint64');
        measLength{k} = fread(fid,1,'uint64'); 
        fseek(fid, 152 - 16, 'cof');
    end
else
    % in VB versions, the first 4 bytes indicate the beginning of the
    % raw data part of the file
    version  = 'vb';
    disp('Software version: VB (!?)');
    measOffset{1} = 0;
    measLength{1} = fileSize;
    NScans     = 1; % VB does not support multiple scans in one file
end

%SRY read data correction factors
% do this for all VB datasets, so that the factors are available later
% in the image_obj if the user chooses to set the correction flag
if (strcmp(version, 'vb')) % not implemented/tested for vd, yet
    datStart = measOffset{1} + firstInt;
    frewind(fid);
    while ( (ftell(fid) < datStart) && ~exist('rawfactors', 'var'))
        line = fgetl(fid);
        %find the section of the protocol
        %note: the factors are also available in <ParamArray."CoilSelects">
        %along with element name and FFT scale, but the parsing is
        %significantly more difficult
        if (~isempty(strfind(line, '<ParamArray."axRawDataCorrectionFactor">')))
            while (ftell(fid) < datStart)
                line = fgetl(fid);
                %find the line with correction factors
                %the factors are on the first line that begins with this
                %pattern
                if (~isempty(strfind(line, '{ {  { ')))
                    line = strrep(line, '}  { } ', '0.0');
                    line = strrep(line, '{', '');
                    line = strrep(line, '}', '');
                    rawfactors = textscan(line, '%f');
                    rawfactors = rawfactors{1}; %textscan returns a 1x1 cell array
                    % this does not work in this location because the MDHs
                    % have not been parsed yet
                    %                    if (length(rawfactors) ~= 2*max(image_obj.NCha))
                    %                       error('Number of raw factors (%f) does not equal channel count (%d)', length(rawfactors)/2, image_obj.NCha);
                    %                    end;
                    if (mod(length(rawfactors),2) ~= 0)
                        error('Error reading rawfactors');
                    end;
                    %note the transpose, this makes the vector
                    %multiplication during the read easier
                    arg.rawDataCorrectionFactors = rawfactors(1:2:end).' + 1i*rawfactors(2:2:end).';
                    break;
                end
            end
        end
    end
    disp('Read raw data correction factors');
end

% data will be read in two steps (two while loops):
%   1) reading all MDHs to find maximum line no., partition no.,... for
%      ima, ref,... scan
%   2) reading the data
tic;
twix_obj        = cell(1,NScans);

for s=1:NScans
    cPos = measOffset{s};
    fseek(fid,cPos,'bof');
    hdr_len = fread(fid, 1,'uint32');

    % read header and calculate regridding (optional)
    rstraj = [];
    if arg.bReadHeader
        [twix_obj{s}.hdr, rstraj] = read_twix_hdr(fid);
    end

    % declare data objects:
    twix_obj{s}.image         = twix_map_obj(arg,'image',filename,version,rstraj);
    twix_obj{s}.noise         = twix_map_obj(arg,'noise',filename,version);
    twix_obj{s}.phasecor      = twix_map_obj(arg,'phasecor',filename,version,rstraj);
    twix_obj{s}.phasestab     = twix_map_obj(arg,'phasestab',filename,version,rstraj);
    twix_obj{s}.phasestabRef0 = twix_map_obj(arg,'phasestab_ref0',filename,version,rstraj);
    twix_obj{s}.phasestabRef1 = twix_map_obj(arg,'phasestab_ref1',filename,version,rstraj);
    twix_obj{s}.refscan       = twix_map_obj(arg,'refscan',filename,version,rstraj);
    twix_obj{s}.refscanPC     = twix_map_obj(arg,'refscan_phasecor',filename,version,rstraj);
    twix_obj{s}.refscanPS     = twix_map_obj(arg,'refscan_phasestab',filename,version,rstraj);
    twix_obj{s}.refscanPSRef0 = twix_map_obj(arg,'refscan_phasestab_ref0',filename,version,rstraj);
    twix_obj{s}.refscanPSRef1 = twix_map_obj(arg,'refscan_phasestab_ref1',filename,version,rstraj);
    twix_obj{s}.RTfeedback    = twix_map_obj(arg,'rtfeedback',filename,version,rstraj);
    twix_obj{s}.vop           = twix_map_obj(arg,'vop',filename,version); % tx-array rf pulses
    
    % print reader version information
    if s==1
        fprintf('Reader version: %d', twix_obj{s}.image.readerVersion);
        isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
        if isOctave
            date_str = ctime(twix_obj{s}.image.readerVersion);
            date_str = date_str(1:end-1);
        else
            date_str = datetime(twix_obj{s}.image.readerVersion, 'ConvertFrom','posixtime');
        end
        fprintf(' (UTC: %s)\n', datestr(date_str));  %Jamie Near - pass date_str as string for 
                                                     %compatibility with older versions of MATLAB.                                                  
    end
    
    % jump to first mdh
    cPos = cPos + hdr_len;
    fseek( fid, cPos, 'bof' );

    % find all mdhs and save them in binary form, first:
    fprintf('Scan %d/%d, read all mdhs:\n', s, NScans )

    [mdh_blob, filePos, isEOF] = loop_mdh_read( fid, version, NScans, s, measOffset{s}, measLength{s});  % uint8; size: [ byteMDH  Nmeas ]

    cPos = filePos( end );
    filePos( end ) = [];

    % get mdhs and masks for each scan, no matter if noise, image, RTfeedback etc:
    [mdh, mask] = evalMDH( mdh_blob, version ); % this is quasi-instant (< 1s) :-)    
    clear mdh_blob

    % Assign mdhs to their respective scans and parse it in the correct twix objects.
    %
    if arg.bReadImaScan
        clear tmpMdh
        isCurrScan = logical( mask.MDH_IMASCAN );
        for f = fieldnames( mdh ).'
            tmpMdh.(f{1}) = mdh.(f{1})( isCurrScan, : );
        end
        twix_obj{s}.image.readMDH( tmpMdh, filePos(isCurrScan) );
    end
    if arg.bReadNoiseScan
        clear tmpMdh
        isCurrScan = logical( mask.MDH_NOISEADJSCAN );
        for f = fieldnames( mdh ).'
            tmpMdh.(f{1}) = mdh.(f{1})( isCurrScan, : );
        end
        twix_obj{s}.noise.readMDH( tmpMdh, filePos(isCurrScan) );
    end
    if arg.bReadRefScan
        clear tmpMdh
        isCurrScan =    ( mask.MDH_PATREFSCAN | mask.MDH_PATREFANDIMASCAN )...
                     & ~( mask.MDH_PHASCOR | mask.MDH_PHASESTABSCAN | ...
                          mask.MDH_REFPHASESTABSCAN | ...
                          mask.MDH_RTFEEDBACK | mask.MDH_HPFEEDBACK);
        for f = fieldnames( mdh ).'
            tmpMdh.(f{1}) = mdh.(f{1})( isCurrScan, : );
        end
        twix_obj{s}.refscan.readMDH( tmpMdh, filePos(isCurrScan) );
    end
    if arg.bReadRTfeedback
        clear tmpMdh
        isCurrScan = ( mask.MDH_RTFEEDBACK | mask.MDH_HPFEEDBACK ) & ~mask.MDH_VOP;
        for f = fieldnames( mdh ).'
            tmpMdh.(f{1}) = mdh.(f{1})( isCurrScan, : );
        end
        twix_obj{s}.RTfeedback.readMDH( tmpMdh, filePos(isCurrScan) );

        clear tmpMdh
        isCurrScan = ( mask.MDH_RTFEEDBACK & mask.MDH_VOP );
        for f = fieldnames( mdh ).'
            tmpMdh.(f{1}) = mdh.(f{1})( isCurrScan, : );
        end
        twix_obj{s}.vop.readMDH( tmpMdh, filePos(isCurrScan) );
    end
    if arg.bReadPCScan
        % logic really correct?
        
        clear tmpMdh
        isCurrScan = mask.MDH_PHASCOR & ( ~mask.MDH_PATREFSCAN | mask.MDH_PATREFANDIMASCAN );
        for f = fieldnames( mdh ).'
            tmpMdh.(f{1}) = mdh.(f{1})( isCurrScan, : );
        end
        twix_obj{s}.phasecor.readMDH( tmpMdh, filePos(isCurrScan) );
        
        clear tmpMdh
        isCurrScan = mask.MDH_PHASCOR & (  mask.MDH_PATREFSCAN | mask.MDH_PATREFANDIMASCAN );
        for f = fieldnames( mdh ).'
            tmpMdh.(f{1}) = mdh.(f{1})( isCurrScan, : );
        end
        twix_obj{s}.refscanPC.readMDH( tmpMdh, filePos(isCurrScan) );
    end
    if arg.bReadPhaseStab
        clear tmpMdh
        isCurrScan =    ( mask.MDH_PHASESTABSCAN & ~mask.MDH_REFPHASESTABSCAN ) ...
                      & (~mask.MDH_PATREFSCAN    |  mask.MDH_PATREFANDIMASCAN );
        for f = fieldnames( mdh ).'
            tmpMdh.(f{1}) = mdh.(f{1})( isCurrScan, : );
        end
        twix_obj{s}.phasestab.readMDH( tmpMdh, filePos(isCurrScan) );

        clear tmpMdh
        isCurrScan =    ( mask.MDH_PHASESTABSCAN & ~mask.MDH_REFPHASESTABSCAN ) ...
                      & ( mask.MDH_PATREFSCAN    |  mask.MDH_PATREFANDIMASCAN );
        for f = fieldnames( mdh ).'
            tmpMdh.(f{1}) = mdh.(f{1})( isCurrScan, : );
        end
        twix_obj{s}.refscanPS.readMDH( tmpMdh, filePos(isCurrScan) );

        clear tmpMdh
        isCurrScan =    ( mask.MDH_REFPHASESTABSCAN & ~mask.MDH_PHASESTABSCAN ) ...
                      & (~mask.MDH_PATREFSCAN   |   mask.MDH_PATREFANDIMASCAN );
        for f = fieldnames( mdh ).'
            tmpMdh.(f{1}) = mdh.(f{1})( isCurrScan, : );
        end
        twix_obj{s}.phasestabRef0.readMDH( tmpMdh, filePos(isCurrScan) );

        clear tmpMdh
        isCurrScan =    ( mask.MDH_REFPHASESTABSCAN & ~mask.MDH_PHASESTABSCAN ) ...
                      & ( mask.MDH_PATREFSCAN   |   mask.MDH_PATREFANDIMASCAN );
        for f = fieldnames( mdh ).'
            tmpMdh.(f{1}) = mdh.(f{1})( isCurrScan, : );
        end
        twix_obj{s}.refscanPSRef0.readMDH( tmpMdh, filePos(isCurrScan) );
        
        clear tmpMdh
        isCurrScan =    ( mask.MDH_REFPHASESTABSCAN & mask.MDH_PHASESTABSCAN ) ...
                      & (~mask.MDH_PATREFSCAN   |   mask.MDH_PATREFANDIMASCAN );
        for f = fieldnames( mdh ).'
            tmpMdh.(f{1}) = mdh.(f{1})( isCurrScan, : );
        end
        twix_obj{s}.phasestabRef1.readMDH( tmpMdh, filePos(isCurrScan) );

        clear tmpMdh
        isCurrScan =    ( mask.MDH_REFPHASESTABSCAN & mask.MDH_PHASESTABSCAN ) ...
                      & ( mask.MDH_PATREFSCAN   |   mask.MDH_PATREFANDIMASCAN );
        for f = fieldnames( mdh ).'
            tmpMdh.(f{1}) = mdh.(f{1})( isCurrScan, : );
        end
        twix_obj{s}.refscanPSRef1.readMDH( tmpMdh, filePos(isCurrScan) );
    end
    clear  mdh  tmpMdh  filePos  isCurrScan

    for scan = { 'image', 'noise', 'phasecor', 'phasestab', ...
                 'phasestabRef0', 'phasestabRef1', 'refscan', ...
                 'refscanPC', 'refscanPS', 'refscanPSRef0', ...
                 'refscanPSRef1', 'RTfeedback', 'vop' }
        f = scan{1};

        % remove unused fields
        if twix_obj{s}.(f).NAcq == 0
            twix_obj{s} = rmfield(twix_obj{s}, f );
        else
            if isEOF
                % recover from read error
                twix_obj{s}.(f).tryAndFixLastMdh();
            else
                twix_obj{s}.(f).clean();
            end
        end
    end

end % NScans loop

if NScans == 1
    twix_obj = twix_obj{1};
end

end % of mapVBVD()




function [mdh_blob, filePos, isEOF] = loop_mdh_read( fid, version, Nscans, scan, measOffset, measLength)
% Goal of this function is to gather all mdhs in the dat file and store them
% in binary form, first. This enables us to evaluate and parse the stuff in
% a MATLAB-friendly (vectorized) way. We also yield a clear separation between
% a lengthy loop and other expressions that are evaluated very few times.
%
% The main challenge is that we never know a priori, where the next mdh is
% and how many there are. So we have to actually evaluate some mdh fields to
% find the next one.
%
% All slow things of the parsing step are found in the while loop.
% => It is the (only) place where micro-optimizations are worthwhile.
%
% The current state is that we are close to sequential disk I/O times.
% More fancy improvements may be possible by using workers through parfeval()
% or threads using a java class (probably faster + no toolbox):
% http://undocumentedmatlab.com/blog/explicit-multi-threading-in-matlab-part1

    switch version
        case 'vb'
            isVD    = false;
            byteMDH = 128;
        case 'vd'
            isVD    = true;
            byteMDH = 184;
            szScanHeader    = 192; % [bytes]
            szChannelHeader =  32; % [bytes]
        otherwise
            % arbitrary assumptions:
            isVD    = false;
            byteMDH = 128;
            warning( [mfilename() ':UnknownVer'], 'Software version "%s" is not supported.', version );
    end

    cPos          = ftell(fid);
    n_acq         = 0;
    allocSize     = 4096;
    ulDMALength   = byteMDH;
    isEOF         = false;
    last_progress = 0;

    mdh_blob = zeros( byteMDH, 0, 'uint8' );
    szBlob   = size( mdh_blob, 2 );
    filePos  = zeros(0, 1, class(cPos));  % avoid bug in Matlab 2013b: https://scivision.co/matlab-fseek-bug-with-uint64-offset/

    fseek(fid,cPos,'bof');

    % ======================================
    %   constants and conditional variables
    % ======================================
        bit_0 = uint8(2^0);
        bit_5 = uint8(2^5);
        mdhStart = 1-byteMDH;
        
        u8_000 = zeros( 3, 1, 'uint8'); % for comparison with data_u8(1:3)

        % 20 fill bytes in VD (21:40)
        evIdx   = uint8(    21  + 20*isVD); % 1st byte of evalInfoMask
        dmaIdx  = uint8((29:32) + 20*isVD); % to correct DMA length using NCol and NCha
        if isVD
            dmaOff  = szScanHeader;
            dmaSkip = szChannelHeader;
        else
            dmaOff  = 0;
            dmaSkip = byteMDH;
        end
    % ======================================

    isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    if isOctave % octave does not support a cancel button
        h = waitbar(0,'','Name', sprintf('Reading Scan ID %d/%d', scan, Nscans));
    else
        h = waitbar(0,'','Name', sprintf('Reading Scan ID %d/%d', scan, Nscans),...
                'CreateCancelBtn',...
                'setappdata(gcbf,''canceling'',1)');
        setappdata(h,'canceling',0)
    end

    t0 = tic;
    while true
        % Read mdh as binary (uint8) and evaluate as little as possible to know...
        %   ... where the next mdh is (ulDMALength / ushSamplesInScan & ushUsedChannels)
        %   ... whether it is only for sync (MDH_SYNCDATA)
        %   ... whether it is the last one (MDH_ACQEND)
        % evalMDH() contains the correct and readable code for all mdh entries.
         
        try
            % read everything and cut out the mdh
            data_u8 = fread( fid, ulDMALength, 'uint8=>uint8' );
            data_u8 = data_u8( mdhStart+end :  end );
        catch exc
            warning( [mfilename() ':UnxpctdEOF'],  ...
                      [ '\nAn unexpected read error occurred at this byte offset: %d (%g GiB)\n'...
                        'Will stop reading now.\n'                                             ...
                        '=== MATLABs error message ================\n'                         ...
                        exc.message                                                            ...
                        '\n=== end of error =========================\n'                       ...
                       ], cPos, cPos/1024^3 )
            isEOF = true;
            break
        end

        if ~isOctave && getappdata(h,'canceling') 
            break;
        end
        
        bitMask = data_u8(evIdx);   % the initial 8 bit from evalInfoMask are enough

        if   isequal( data_u8(1:3), u8_000 )    ... % probably ulDMALength == 0
          || bitand(bitMask, bit_0);                % MDH_ACQEND

            % ok, look closer if really all *4* bytes are 0:
            data_u8(4)= bitget( data_u8(4),1);  % ubit24: keep only 1 bit from the 4th byte
            ulDMALength = double( typecast( data_u8(1:4), 'uint32' ) );

            if ulDMALength == 0 || bitand(bitMask, bit_0)
                cPos = cPos + ulDMALength;
                % jump to next full 512 bytes
                if mod(cPos,512)
                    cPos = cPos + 512 - mod(cPos,512);
                end
                break;
            end
        end
        if bitand(bitMask, bit_5);  % MDH_SYNCDATA
            data_u8(4)= bitget( data_u8(4),1);  % ubit24: keep only 1 bit from the 4th byte
            ulDMALength = double( typecast( data_u8(1:4), 'uint32' ) );
            cPos = cPos + ulDMALength;
            continue
        end

        % pehses: the pack bit indicates that multiple ADC are packed into one
        % DMA, often in EPI scans (controlled by fRTSetReadoutPackaging in IDEA)
        % since this code assumes one adc (x NCha) per DMA, we have to correct
        % the "DMA length"
        %     if mdh.ulPackBit
        % it seems that the packbit is not always set correctly
        NCol_NCha = double( typecast( data_u8(dmaIdx), 'uint16' ) );  % [ushSamplesInScan  ushUsedChannels]
        ulDMALength = dmaOff + (8*NCol_NCha(1) + dmaSkip) * NCol_NCha(2);

        n_acq = n_acq + 1;

        % grow arrays in batches
        if n_acq > szBlob
            mdh_blob( :, end + allocSize ) = 0;
            filePos( end + allocSize ) = 0;
            szBlob = size( mdh_blob, 2 );
        end
        mdh_blob(:,n_acq) = data_u8;
        filePos( n_acq )  = cPos;

        progress = (cPos-measOffset)/measLength;
        
        if progress > last_progress  + 0.01
            last_progress = progress;
            elapsed_time  = toc(t0);
            time_left     = elapsed_time * (1/progress-1);
            progress_str  = sprintf('%3.0f %% read in %4.0f s;\nestimated time left: %4.0f s', round(100*progress), elapsed_time, time_left);
            waitbar(progress, h, progress_str);
        end

        cPos = cPos + ulDMALength;
    end % while true
    delete(h);
    
    if isEOF
        n_acq = n_acq-1;    % ignore the last attempt
    end

    filePos( n_acq+1 ) = cPos;  % save pointer to the next scan

    % discard overallocation:
    mdh_blob = mdh_blob(:,1:n_acq);
    filePos  = reshape( filePos(1:n_acq+1), 1, [] ); % row vector

    fprintf('%8.1f MB read in %4.0f s\n', measLength/1024^2, round(toc(t0)));

end % of loop_mdh_read()



function [mdh,mask] = evalMDH( mdh_blob, version )
% see pkg/MrServers/MrMeasSrv/SeqIF/MDH/mdh.h
% and pkg/MrServers/MrMeasSrv/SeqIF/MDH/MdhProxy.h

if ~isa( mdh_blob, 'uint8' )
    error([mfilename() ':NoInt8'], 'mdh data must be a uint8 array!')
end

if version(end) == 'd'
    isVD = true;
    mdh_blob = mdh_blob([1:20 41:end], :);  % remove 20 unnecessary bytes
else
    isVD = false;
end

Nmeas   = size( mdh_blob, 2 );

mdh.ulPackBit   = bitget( mdh_blob(4,:), 2).';
mdh.ulPCI_rx    = bitset(bitset(mdh_blob(4,:), 7, 0), 8, 0).'; % keep 6 relevant bits
mdh_blob(4,:)   = bitget( mdh_blob(4,:),1);  % ubit24: keep only 1 bit from the 4th byte

% unfortunately, typecast works on vectors, only
data_uint32     = typecast( reshape(mdh_blob(1:76,:),  [],1), 'uint32' );
data_uint16     = typecast( reshape(mdh_blob(29:end,:),[],1), 'uint16' );
data_single     = typecast( reshape(mdh_blob(69:end,:),[],1), 'single' );

data_uint32 = reshape( data_uint32, [], Nmeas ).';
data_uint16 = reshape( data_uint16, [], Nmeas ).';
data_single = reshape( data_single, [], Nmeas ).';
                                                        %  byte pos.
%mdh.ulDMALength               = data_uint32(:,1);      %   1 :   4
mdh.lMeasUID                   = data_uint32(:,2);      %   5 :   8
mdh.ulScanCounter              = data_uint32(:,3);      %   9 :  12
mdh.ulTimeStamp                = data_uint32(:,4);      %  13 :  16
mdh.ulPMUTimeStamp             = data_uint32(:,5);      %  17 :  20
mdh.aulEvalInfoMask            = data_uint32(:,6:7);    %  21 :  28
mdh.ushSamplesInScan           = data_uint16(:,1);      %  29 :  30
mdh.ushUsedChannels            = data_uint16(:,2);      %  31 :  32
mdh.sLC                        = data_uint16(:,3:16);   %  33 :  60
mdh.sCutOff                    = data_uint16(:,17:18);  %  61 :  64
mdh.ushKSpaceCentreColumn      = data_uint16(:,19);     %  66 :  66
mdh.ushCoilSelect              = data_uint16(:,20);     %  67 :  68
mdh.fReadOutOffcentre          = data_single(:, 1);     %  69 :  72
mdh.ulTimeSinceLastRF          = data_uint32(:,19);     %  73 :  76
mdh.ushKSpaceCentreLineNo      = data_uint16(:,25);     %  77 :  78
mdh.ushKSpaceCentrePartitionNo = data_uint16(:,26);     %  79 :  80

if isVD
    mdh.SlicePos                    = data_single(:, 4:10); %  81 : 108
    mdh.aushIceProgramPara          = data_uint16(:,41:64); % 109 : 156
    mdh.aushFreePara                = data_uint16(:,65:68); % 157 : 164
else
    mdh.aushIceProgramPara          = data_uint16(:,27:30); %  81 :  88
    mdh.aushFreePara                = data_uint16(:,31:34); %  89 :  96
    mdh.SlicePos                    = data_single(:, 8:14); %  97 : 124
end

% inlining of evalInfoMask
evalInfoMask1 = mdh.aulEvalInfoMask(:,1);
mask.MDH_ACQEND            = min(bitand(evalInfoMask1, 2^0), 1);
mask.MDH_RTFEEDBACK        = min(bitand(evalInfoMask1, 2^1), 1);
mask.MDH_HPFEEDBACK        = min(bitand(evalInfoMask1, 2^2), 1);
mask.MDH_SYNCDATA          = min(bitand(evalInfoMask1, 2^5), 1);
mask.MDH_RAWDATACORRECTION = min(bitand(evalInfoMask1, 2^10),1);
mask.MDH_REFPHASESTABSCAN  = min(bitand(evalInfoMask1, 2^14),1);
mask.MDH_PHASESTABSCAN     = min(bitand(evalInfoMask1, 2^15),1);
mask.MDH_SIGNREV           = min(bitand(evalInfoMask1, 2^17),1);
mask.MDH_PHASCOR           = min(bitand(evalInfoMask1, 2^21),1);
mask.MDH_PATREFSCAN        = min(bitand(evalInfoMask1, 2^22),1);
mask.MDH_PATREFANDIMASCAN  = min(bitand(evalInfoMask1, 2^23),1);
mask.MDH_REFLECT           = min(bitand(evalInfoMask1, 2^24),1);
mask.MDH_NOISEADJSCAN      = min(bitand(evalInfoMask1, 2^25),1);
mask.MDH_VOP               = min(bitand(mdh.aulEvalInfoMask(2), 2^(53-32)),1); % was 0 in VD
mask.MDH_IMASCAN           = ones( Nmeas, 1, 'uint32' );

noImaScan = ( mask.MDH_ACQEND           | mask.MDH_RTFEEDBACK   | mask.MDH_HPFEEDBACK    ...
            | mask.MDH_PHASCOR          | mask.MDH_NOISEADJSCAN | mask.MDH_PHASESTABSCAN ...
            | mask.MDH_REFPHASESTABSCAN | mask.MDH_SYNCDATA                              ... 
            | (mask.MDH_PATREFSCAN & ~mask.MDH_PATREFANDIMASCAN) );

mask.MDH_IMASCAN( noImaScan ) = 0;

end % of evalMDH()
