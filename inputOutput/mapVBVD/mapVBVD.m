function twix_obj = mapVBVD(filename,varargin)

%  Reads Siemens raw .dat file from VB/VD MRI raw data.
%
%  Requires twix_map_obj.m
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
% 
% Input:
%
% filename or simply meas. id, e.g. mapVBVD(122) (if file is in same path)
% optional arguments (see below)
%
% Output: twix_obj structure with elements (if available):
%       .image:       object for image scan
%       .noise:       object for noise scan
%       .phasecor:    object for phase correction scan
%       .refscan:     object for reference scan
%       .refscanPC:   object for phase correction scan for reference data
%       .RTfeedback:  object for realtime feedback data
%       .phasestab:   object for phase stabilization scan
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
% noWeightedAverage: disables proper weighting of averages (i.e. data is 
%                    summed up instead of averaged)
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
        if regexp(filelist(k).name,['^meas_MID0*' num2str(measID) '.*\.dat'])==1
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


% add absolute path, when no path is given
[pathstr, name, ext] = fileparts(filename);

if isempty(pathstr)
    pathstr  = pwd;
    filename = fullfile(pathstr, [name ext]);
end

%%%%% Parse varargin %%%%%

    % Definition of default parameters
    arg.removeOS          = false;
    arg.doAverage         = false;
    arg.averageReps       = false;
    arg.averageSets       = false;
    arg.ignoreSeg         = false;
    arg.noWeightedAverage = false;
    arg.bReadImaScan      = true;
    arg.bReadNoiseScan    = true;
    arg.bReadPCScan       = true;
    arg.bReadRefScan      = true;
    arg.bReadRefPCScan    = true;
    arg.bReadRTfeedback   = true;
    arg.bReadPhaseStab    = true;
    arg.doRawDataCorrect  = false; %SRY
    
    k=1;
    while k <= numel(varargin)
        
        if ~ischar(varargin{k})
            error('string expected');
        end
        
        switch lower(varargin{k})
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
            case {'noweightedaverage','noweightedave','noweighting'}
                if numel(varargin) > k && ~ischar(varargin{k+1})
                    arg.noWeightedAverage = logical(varargin{k+1});
                    k = k+2;
                else
                    arg.noWeightedAverage = true;
                    k = k+1;
                end
            case {'rawdatacorrect','dorawdatacorrect'}
                %SRY: handle raw data correct arguments
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
    fid = fopen(filename,'r','l','US-ASCII');
    fseek(fid,0,'eof');
    fileSize = ftell(fid);
    
    % start of actual measurment data (sans header)
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
        % measOffset: points to beginning of header, usually at 10240 bytes
        measOffset = fread(fid,1,'uint64');
        measLength = fread(fid,1,'uint64');
        fseek(fid,measOffset,'bof');
        hdrLength  = fread(fid,1,'uint32');
        datStart   = measOffset + hdrLength;
    else
        % in VB versions, the first 4 bytes indicate the beginning of the
        % raw data part of the file
        version  = 'vb';
        disp('Software version: VB (!?)');
        datStart = firstInt;
        NScans   = 1; % VB does not support multiple scans in one file
    end
    
    %SRY read data correction factors
    % do this for all VB datasets, so that the factors are available later
    % in the image_obj if the user chooses to set the correction flag
    if (strcmp(version, 'vb')) % not implemented/tested for vd, yet
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
    percentFinished = 0;
    cPos            = datStart;
    twix_obj = cell(1,NScans);
    for s=1:NScans
        % declare data objects:
        twix_obj{s}.image      = twix_map_obj(arg,'image',filename,version);
        twix_obj{s}.noise      = twix_map_obj(arg,'noise',filename,version);
        twix_obj{s}.phasecor   = twix_map_obj(arg,'phasecor',filename,version);
        twix_obj{s}.refscan    = twix_map_obj(arg,'refscan',filename,version);
        twix_obj{s}.refscanPC  = twix_map_obj(arg,'refscan_phasecor',filename,version);
        twix_obj{s}.RTfeedback = twix_map_obj(arg,'rtfeedback',filename,version);    
        twix_obj{s}.phasestab  = twix_map_obj(arg,'phasestab',filename,version);

        mask.MDH_ACQEND = 0;
        frewind(fid);
        while cPos+128<fileSize % fail-safe; in case we miss MDH_ACQEND

            switch version
                case 'vb'
                    [mdh mask] = evalMDHvb(fid,cPos);
                case 'vd'
                    [mdh mask] = evalMDHvd(fid,cPos);
                otherwise
                    disp('error: only vb/vd software versions supported');                    
            end

            if mask.MDH_ACQEND || mdh.ulDMALength==0
                if s<NScans
                    cPos = cPos + mdh.ulDMALength;
                    % jump to next full 512 bytes
                    cPos = cPos + 512 - mod(cPos,512);
                    fseek(fid,cPos,'bof');
                    hdrLength  = fread(fid,1,'uint32');
                    cPos = cPos + hdrLength;
                end
                break;
            end

            if mask.MDH_SYNCDATA
                % skip SYNCDATA
                cPos = cPos + mdh.ulDMALength;
                continue;
            end

            if (mask.MDH_IMASCAN && arg.bReadImaScan)
                twix_obj{s}.image.readMDH(mdh,cPos);
            end

            if (mask.MDH_NOISEADJSCAN && arg.bReadNoiseScan)
                twix_obj{s}.noise.readMDH(mdh,cPos);
            end

            if (and(mask.MDH_PHASCOR,~mask.MDH_PATREFSCAN) && arg.bReadPCScan)
                twix_obj{s}.phasecor.readMDH(mdh,cPos);
            end

            if (and(~mask.MDH_PHASCOR,(mask.MDH_PATREFSCAN || mask.MDH_PATREFANDIMASCAN)) && arg.bReadRefScan)
                twix_obj{s}.refscan.readMDH(mdh,cPos);
            end

            if (and(mask.MDH_PATREFSCAN,mask.MDH_PHASCOR) && arg.bReadRefPCScan)
                twix_obj{s}.refscanPC.readMDH(mdh,cPos);
            end

            if ((mask.MDH_RTFEEDBACK || mask.MDH_HPFEEDBACK) && arg.bReadRTfeedback)
                twix_obj{s}.RTfeedback.readMDH(mdh,cPos);
            end

            if ((mask.MDH_PHASESTABSCAN || mask.MDH_REFPHASESTABSCAN) && arg.bReadPhaseStab)
                twix_obj{s}.phasestab.readMDH(mdh,cPos);
            end

            % jump to mdh of next scan
            cPos = cPos + mdh.ulDMALength;

            if (cPos/fileSize*100 > percentFinished + 1)
                percentFinished = floor(cPos/fileSize*100);
                elapsed_time  = toc;
                time_left     = (fileSize/cPos-1) * elapsed_time;

                if ~exist('progress_str','var')
                    prevLength = 0;
                else
                    prevLength = numel(progress_str);
                end

                progress_str = sprintf('%3.0f %% parsed in %4.0f s; estimated time left: %4.0f s \n',...
                percentFinished,elapsed_time, time_left);

                fprintf([repmat('\b',1,prevLength) '%s'],progress_str);
            end

        end % while         

        if twix_obj{s}.image.NAcq == 0
            twix_obj{s} = rmfield(twix_obj{s},'image');
        else
            twix_obj{s}.image.clean();
        end

        if twix_obj{s}.noise.NAcq == 0
            twix_obj{s} = rmfield(twix_obj{s},'noise');
        else
            twix_obj{s}.noise.clean();
        end

        if twix_obj{s}.phasecor.NAcq == 0
            twix_obj{s} = rmfield(twix_obj{s},'phasecor');
        else
            twix_obj{s}.phasecor.clean();
        end

        if twix_obj{s}.refscan.NAcq == 0
            twix_obj{s} = rmfield(twix_obj{s},'refscan');
        else
            twix_obj{s}.refscan.clean();
        end

        if twix_obj{s}.refscanPC.NAcq == 0
            twix_obj{s} = rmfield(twix_obj{s},'refscanPC');
        else
            twix_obj{s}.refscanPC.clean();
        end

        if twix_obj{s}.RTfeedback.NAcq == 0
            twix_obj{s} = rmfield(twix_obj{s},'RTfeedback');
        else
            twix_obj{s}.RTfeedback.clean();
        end

        if twix_obj{s}.phasestab.NAcq == 0
            twix_obj{s} = rmfield(twix_obj{s},'phasestab');
        else
            twix_obj{s}.phasestab.clean();
        end
    end % NScans loop
   
    if NScans == 1
        twix_obj = twix_obj{1};
    end
    
    elapsed_time = toc;
    progress_str = sprintf('100 %% parsed in %4.0f s; estimated time left:    0 s \n', elapsed_time);
    fprintf([repmat('\b',1,prevLength) '%s'],progress_str);
end


function [mdh mask] = evalMDHvb(fid,cPos)
    % see pkg/MrServers/MrMeasSrv/SeqIF/MDH/mdh.h
    % and pkg/MrServers/MrMeasSrv/SeqIF/MDH/MdhProxy.h
    
    % no difference between 'scan' and 'channel' header in VB
    szMDH = 128; % [bytes]

    % inlining of readMDH
    fseek(fid,cPos,'bof');
    mdh.ulDMALength                = fread(fid,      1, 'ubit25');
    mdh.ulPackBit                  = fread(fid,      1, 'ubit1');
    mdh.ulPCI_rx                   = fread(fid,      1, 'ubit6');
    
    fseek(fid,cPos,'bof');
    mdh.ulFlagsAndDMALength        = fread(fid,      1, 'uint32');
    fseek(fid,cPos+20,'bof');
    mdh.aulEvalInfoMask            = fread(fid,  [1 2], 'uint32');
    dummy                          = fread(fid,      2, 'uint16');
    mdh.ushSamplesInScan           = dummy(1);
    mdh.ushUsedChannels            = dummy(2);
    mdh.sLC                        = fread(fid, [1 14], 'ushort');  
    dummy                          = fread(fid,     18, 'uint16');
%     mdh.sCutOff                    = dummy(1:2);
    mdh.ushKSpaceCentreColumn      = dummy(3);
%     mdh.ushCoilSelect              = dummy(4);
    mdh.ushKSpaceCentreLineNo      = dummy(9);
    mdh.ushKSpaceCentrePartitionNo = dummy(10);
    mdh.aushIceProgramPara         = dummy(11:14);
    mdh.aushFreePara               = dummy(15:18);
    mdh.SlicePos                   = fread(fid,      7, 'float');
    
    % inlining of evalInfoMask
    mask.MDH_ACQEND             = min(bitand(mdh.aulEvalInfoMask(1), 2^0),1);
    mask.MDH_RTFEEDBACK         = min(bitand(mdh.aulEvalInfoMask(1), 2^1),1);
    mask.MDH_HPFEEDBACK         = min(bitand(mdh.aulEvalInfoMask(1), 2^2),1);
    mask.MDH_SYNCDATA           = min(bitand(mdh.aulEvalInfoMask(1), 2^5), 1);
    mask.MDH_RAWDATACORRECTION  = min(bitand(mdh.aulEvalInfoMask(1), 2^10),1);
    mask.MDH_REFPHASESTABSCAN   = min(bitand(mdh.aulEvalInfoMask(1), 2^14),1);
    mask.MDH_PHASESTABSCAN      = min(bitand(mdh.aulEvalInfoMask(1), 2^15),1);
    mask.MDH_SIGNREV            = min(bitand(mdh.aulEvalInfoMask(1), 2^17),1);
    mask.MDH_PHASCOR            = min(bitand(mdh.aulEvalInfoMask(1), 2^21),1);
    mask.MDH_PATREFSCAN         = min(bitand(mdh.aulEvalInfoMask(1), 2^22),1);
    mask.MDH_PATREFANDIMASCAN   = min(bitand(mdh.aulEvalInfoMask(1), 2^23),1);
    mask.MDH_REFLECT            = min(bitand(mdh.aulEvalInfoMask(1), 2^24),1);
    mask.MDH_NOISEADJSCAN       = min(bitand(mdh.aulEvalInfoMask(1), 2^25),1);
    
    mask.MDH_IMASCAN            = 1;
    if (mask.MDH_ACQEND || mask.MDH_RTFEEDBACK || mask.MDH_HPFEEDBACK...
                        || mask.MDH_PHASCOR    || mask.MDH_NOISEADJSCAN...
                        || mask.MDH_SYNCDATA)
        mask.MDH_IMASCAN = 0; 
    end
    
    % otherwise the PATREFSCAN may be overwritten
    if mask.MDH_PHASESTABSCAN || mask.MDH_REFPHASESTABSCAN
        mask.MDH_PATREFSCAN = 0;
        mask.MDH_PATREFANDIMASCAN = 0;
        mask.MDH_IMASCAN = 0; 
    end
    
    if ( mask.MDH_PATREFSCAN && ~mask.MDH_PATREFANDIMASCAN )
        mask.MDH_IMASCAN = 0;
    end
    
    % pehses: the pack bit indicates that multiple ADC are packed into one 
    % DMA, often in EPI scans (controlled by fRTSetReadoutPackaging in IDEA)
    % since this code assumes one adc (x NCha) per DMA, we have to correct 
    % the "DMA length"
%     if mdh.ulPackBit
    % it seems that the packbit is not always set correctly
    if ~mask.MDH_SYNCDATA && ~mask.MDH_ACQEND
        mdh.ulDMALength = (2*4*mdh.ushSamplesInScan + szMDH) * mdh.ushUsedChannels;
    end

end
    

function [mdh mask] = evalMDHvd(fid,cPos)
    % see pkg/MrServers/MrMeasSrv/SeqIF/MDH/mdh.h
    % and pkg/MrServers/MrMeasSrv/SeqIF/MDH/MdhProxy.h
    
    % we need to differentiate between 'scan header' and 'channel header'
    % since these are used in VD versions:
    szScanHeader    = 192; % [bytes]
    szChannelHeader = 32;  % [bytes]

    % inlining of readScanHeader
    fseek(fid,cPos,'bof');
    mdh.ulDMALength                = fread(fid,      1, 'ubit25');
    mdh.ulPackBit                  = fread(fid,      1, 'ubit1');
    mdh.ulPCI_rx                   = fread(fid,      1, 'ubit6');
    fseek(fid,cPos+40,'bof');
    mdh.aulEvalInfoMask            = fread(fid,  [1 2], 'uint32');
    dummy                          = fread(fid,      2, 'uint16');
    mdh.ushSamplesInScan           = dummy(1);
    mdh.ushUsedChannels            = dummy(2);
    mdh.sLC                        = fread(fid, [1 14], 'ushort');
    dummy                          = fread(fid,     10, 'uint16');
%     mdh.sCutOff                    = dummy(1:2);
    mdh.ushKSpaceCentreColumn      = dummy(3);
%     mdh.ushCoilSelect              = dummy(4);
    mdh.ushKSpaceCentreLineNo      = dummy(9);
    mdh.ushKSpaceCentrePartitionNo = dummy(10);
    mdh.SlicePos                   = fread(fid,      7, 'float');
    dummy                          = fread(fid,     28, 'uint16');
    mdh.aushIceProgramPara         = dummy(1:24);
    mdh.aushFreePara               = dummy(25:28); % actually aushReservedPara;
                                                   % there's no freePara in VD
  
    % inlining of evalInfoMask
    mask.MDH_ACQEND             = min(bitand(mdh.aulEvalInfoMask(1), 2^0),1);
    mask.MDH_RTFEEDBACK         = min(bitand(mdh.aulEvalInfoMask(1), 2^1),1);
    mask.MDH_HPFEEDBACK         = min(bitand(mdh.aulEvalInfoMask(1), 2^2),1);
    mask.MDH_SYNCDATA           = min(bitand(mdh.aulEvalInfoMask(1), 2^5), 1);
    mask.MDH_RAWDATACORRECTION  = min(bitand(mdh.aulEvalInfoMask(1), 2^10),1);
    mask.MDH_REFPHASESTABSCAN   = min(bitand(mdh.aulEvalInfoMask(1), 2^14),1);
    mask.MDH_PHASESTABSCAN      = min(bitand(mdh.aulEvalInfoMask(1), 2^15),1);
    mask.MDH_SIGNREV            = min(bitand(mdh.aulEvalInfoMask(1), 2^17),1);
    mask.MDH_PHASCOR            = min(bitand(mdh.aulEvalInfoMask(1), 2^21),1);
    mask.MDH_PATREFSCAN         = min(bitand(mdh.aulEvalInfoMask(1), 2^22),1);
    mask.MDH_PATREFANDIMASCAN   = min(bitand(mdh.aulEvalInfoMask(1), 2^23),1);
    mask.MDH_REFLECT            = min(bitand(mdh.aulEvalInfoMask(1), 2^24),1);
    mask.MDH_NOISEADJSCAN       = min(bitand(mdh.aulEvalInfoMask(1), 2^25),1);
    
    mask.MDH_IMASCAN            = 1;
    if (mask.MDH_ACQEND || mask.MDH_RTFEEDBACK || mask.MDH_HPFEEDBACK...
                        || mask.MDH_PHASCOR    || mask.MDH_NOISEADJSCAN...
                        || mask.MDH_SYNCDATA)
        mask.MDH_IMASCAN = 0; 
    end
    
    % otherwise the PATREFSCAN may be overwritten
    if mask.MDH_PHASESTABSCAN || mask.MDH_REFPHASESTABSCAN
        mask.MDH_PATREFSCAN = 0;
        mask.MDH_PATREFANDIMASCAN = 0;
        mask.MDH_IMASCAN = 0; 
    end
    
    if ( mask.MDH_PATREFSCAN && ~mask.MDH_PATREFANDIMASCAN )
        mask.MDH_IMASCAN = 0;
    end
    
    % pehses: the pack bit indicates that multiple ADC are packed into one 
    % DMA, often in EPI scans (controlled by fRTSetReadoutPackaging in IDEA)
    % since this code assumes one adc (x NCha) per DMA, we have to correct 
    % the "DMA length"
%     if mdh.ulPackBit
    % it seems that the packbit is not always set correctly
    if ~mask.MDH_SYNCDATA && ~mask.MDH_ACQEND && mdh.ulDMALength~=0
        mdh.ulDMALength = szScanHeader + (2*4*mdh.ushSamplesInScan + szChannelHeader) * mdh.ushUsedChannels;
    end

end
