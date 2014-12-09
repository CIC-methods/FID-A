classdef twix_map_obj < handle
% class to hold information about raw data from siemens MRI scanners 
% (currently VB and VD software versions are supported and tested).
%
% Author: Philipp Ehses (philipp.ehses@tuebingen.mpg.de), Aug/19/2011
%
%
% Modified by Wolf Blecher (wolf.blecher@tuebingen.mpg.de), Apr/26/2012
% Added reorder index to indicate which lines are reflected
% Added slice position for sorting, Mai/15/2012
%
% Order of many mdh parameters are now stored (including the reflected ADC 
% bit); PE, Jun/14/2012
%
% data is now 'memory mapped' and not read until demanded;
% (see mapVBVD for a description) PE, Aug/02/2012
%
% twix_obj.image.unsorted now returns the data in its acq. order 
% [NCol,NCha,nsamples in acq. order], all average flags don't have an
% influence on the output, but 'flagRemoveOS' still works, PE, Sep/04/13
%
    properties(Dependent=true)
        % flags
        flagRemoveOS     % removes oversampling in read (col) during read operation
        flagDoAverage    % averages over all avg during read operation
        flagAverageReps  % averages over all repetitions
        flagAverageSets  % averages over all sets
    end
    
    properties
        flagNoWeightedAverage % scaling/weighting of averages is disabled
    end
    
    properties(Dependent=true)
        flagIgnoreSeg       % sum over all segments during read operation
        
        flagSkipToFirstLine % skips lines/partitions up to the first 
                            % actually acquired line/partition
                            % (e.g. only the center k-space is acquired in 
                            % refscans, we don't want all the leading zeros
                            % in our data)
                            % this is the default behaviour for everything
                            % but image scans (but can be changed manually)

        flagDoRawDataCorrect     %SRY: apply raw data correction factors during read operation
        RawDataCorrectionFactors %SRY: allow the user to set/get the factors
    end
    
    properties(GetAccess='public', SetAccess='protected')
        % properties:
        filename
        softwareVersion
        dataType
        
        dataSize % this is the current output size, depends on fullSize + some flags
        dataDims
        sqzSize
        sqzDims
        
        NCol  % mdh information
        NCha  % mdh information
        NLin  % mdh information
        NPar  % mdh information
        NSli  % mdh information
        NAve  % mdh information
        NPhs  % mdh information
        NEco  % mdh information
        NRep  % mdh information
        NSet  % mdh information
        NSeg  % mdh information
        NIda  % mdh information
        NIdb  % mdh information
        NIdc  % mdh information
        NIdd  % mdh information
        NIde  % mdh information
        NAcq  % simple counter
        
        % mdh information
        Lin
        Par
        Sli
        Ave
        Phs
        Eco
        Rep
        Set
        Seg
        Ida
        Idb
        Idc
        Idd
        Ide
        
        centerCol
        centerLin
        centerPar
        IsReflected
        IsRawDataCorrect %SRY: storage for MDH flag raw data correct
        
        slicePos
        freeParam
        iceParam
        
        % memory position in file
        memPos
        
        % index that translates simple, linear order of mdh info vectors 
        % to target matrix (of size dataSize)
        ixToTarget % inverted page table (physical to virtual addresses)
        ixToRaw    % page table (virtual to physical addresses)
    end
    
    properties(GetAccess='protected', SetAccess='protected')
        arg  % arguments
        
        allocSize    % determines size of allocation
        currentAlloc % simple counter, keeps track of allocated memory
                
        fullSize % this is the full size of the data set according to the mdhs, i.e. flags
                 % like 'reduceOS' have no influence on it
                 
        freadInfo
        
        skipLin
        skipPar
    end
        
    methods
        % Constructor:
        function this = twix_map_obj(arg,dataType,fname,version)
            
            if ~exist('dataType','var')
                this.dataType = 'image';
            else
                this.dataType = lower(dataType);
            end
            
            this.filename         = fname;
            this.softwareVersion  = version;
            
            this.IsReflected      = logical([]);
            this.IsRawDataCorrect = logical([]); %SRY
            this.NAcq             = 0;
            this.allocSize        = 4096;
            this.currentAlloc     = 0;
                
            if ~isfield(arg,'skipToFirstLine')
                if strcmp(this.dataType,'image')
                    arg.skipToFirstLine = false;
                else
                    arg.skipToFirstLine = true;
                end
            end
                
            if ~exist('arg','var')
                this.arg = [];
            else
                this.arg = arg;
            end
            
            % flags:
            this.flagNoWeightedAverage = arg.noWeightedAverage;
                
            switch this.softwareVersion
                case 'vb'
                    % every channel has its own full mdh
                    this.freadInfo.szScanHeader    =   0; % [bytes]
                    this.freadInfo.szChannelHeader = 128; % [bytes] 
                    this.freadInfo.iceParamSz      =   4;
                case 'vd'
                    if ( this.arg.doRawDataCorrect )
                        error('raw data correction for VD not supported/tested yet');
                    end 
                    this.freadInfo.szScanHeader    = 192; % [bytes]
                    this.freadInfo.szChannelHeader =  32; % [bytes]
                    this.freadInfo.iceParamSz      =  24; % vd version supports up to 24 ice params
                   
                otherwise
                    error('software version not supported');
            end
            
        end
        
            
        function this = readMDH(this,mdh,filePos)

            cLin      = mdh.sLC(1)  + 1;         %%% current line
            cPar      = mdh.sLC(4)  + 1;         %%% current partition  
            cSli      = mdh.sLC(3)  + 1;         %%% current slice
            cAve      = mdh.sLC(2)  + 1;         %%% current scan ('average')
            cPhs      = mdh.sLC(6)  + 1;         %%% current phase cycling step
            cEco      = mdh.sLC(5)  + 1;         %%% current echo no (untested)
            cRep      = mdh.sLC(7)  + 1;         %%% current measurement no
            cSet      = mdh.sLC(8)  + 1;         %%% current set no
            cSeg      = mdh.sLC(9)  + 1;         %%% current segment for future use
            cIda      = mdh.sLC(10) + 1;         %%% ICE dim a
            cIdb      = mdh.sLC(11) + 1;         %%% ICE dim b
            cIdc      = mdh.sLC(12) + 1;         %%% ICE dim c
            cIdd      = mdh.sLC(13) + 1;         %%% ICE dim d
            cIde      = mdh.sLC(14) + 1;         %%% ICE dim e
            
            % subsref overloading makes this.that-calls slow, so we need to
            % avoid them whenever possible
            cAcq      = this.NAcq   + 1;
            this.NAcq = cAcq; 
            
            if cAcq > this.currentAlloc
                % we need to allocate more memory...
                this.currentAlloc = this.currentAlloc + this.allocSize;
                alloc             = zeros(1,this.allocSize  ,'single');
                this.NCol         = cat(2, this.NCol        , alloc);
                this.NCha         = cat(2, this.NCha        , alloc);
                this.Lin          = cat(2, this.Lin         , alloc);
                this.Par          = cat(2, this.Par         , alloc);
                this.Sli          = cat(2, this.Sli         , alloc);
                this.Ave          = cat(2, this.Ave         , alloc);
                this.Phs          = cat(2, this.Phs         , alloc);
                this.Eco          = cat(2, this.Eco         , alloc);
                this.Rep          = cat(2, this.Rep         , alloc);
                this.Set          = cat(2, this.Set         , alloc);
                this.Seg          = cat(2, this.Seg         , alloc);
                this.Ida          = cat(2, this.Ida         , alloc);
                this.Idb          = cat(2, this.Idb         , alloc);
                this.Idc          = cat(2, this.Idc         , alloc);
                this.Idd          = cat(2, this.Idd         , alloc);
                this.Ide          = cat(2, this.Ide         , alloc);
                this.centerCol    = cat(2, this.centerCol   , alloc);
                this.centerLin    = cat(2, this.centerLin   , alloc);
                this.centerPar    = cat(2, this.centerPar   , alloc);
                this.IsReflected  = cat(2, this.IsReflected , false(1,this.allocSize));
                this.IsRawDataCorrect = cat(2, this.IsRawDataCorrect, false(1, this.allocSize)); %SRY
                this.slicePos     = cat(2, this.slicePos    , zeros(7,this.allocSize,'single'));
                this.iceParam     = cat(2, this.iceParam    , zeros(this.freadInfo.iceParamSz, this.allocSize,'single'));
                this.freeParam    = cat(2, this.freeParam   , zeros(4, this.allocSize,'single'));
                this.memPos       = cat(2, this.memPos      , zeros(1,this.allocSize,'double'));
            end
            
            % save mdh information about current line
            this.NCol       (cAcq) = mdh.ushSamplesInScan + 0; % +0 significantly faster??!
            this.NCha       (cAcq) = mdh.ushUsedChannels  + 0;
            this.Lin        (cAcq) = cLin;
            this.Par        (cAcq) = cPar;
            this.Sli        (cAcq) = cSli;
            this.Ave        (cAcq) = cAve;
            this.Phs        (cAcq) = cPhs;
            this.Eco        (cAcq) = cEco;
            this.Rep        (cAcq) = cRep;
            this.Set        (cAcq) = cSet;
            this.Seg        (cAcq) = cSeg;
            this.Ida        (cAcq) = cIda;
            this.Idb        (cAcq) = cIdb;
            this.Idc        (cAcq) = cIdc;
            this.Idd        (cAcq) = cIdd;
            this.Ide        (cAcq) = cIde;
            this.centerCol  (cAcq) = mdh.ushKSpaceCentreColumn + 1;
            this.centerLin  (cAcq) = mdh.ushKSpaceCentreLineNo + 1;
            this.centerPar  (cAcq) = mdh.ushKSpaceCentrePartitionNo + 1;
            this.IsReflected(cAcq) = logical(min(bitand(mdh.aulEvalInfoMask(1),2^24),1));
            this.IsRawDataCorrect(cAcq) = logical(min(bitand(mdh.aulEvalInfoMask(1),2^10),1)); %SRY
            this.slicePos (:,cAcq) = mdh.SlicePos + 0;
            this.iceParam (:,cAcq) = mdh.aushIceProgramPara + 0;
            this.freeParam(:,cAcq) = mdh.aushFreePara + 0;
            % save memory position
            this.memPos     (cAcq) = filePos;
        end
        
        
        function this = clean(this)

            if this.NAcq == 0
                return;
            end
            
            % cut mdh data to actual size (remove over-allocated part)
            this.NCol        = this.NCol       (1:this.NAcq);
            this.NCha        = this.NCha       (1:this.NAcq);
            this.Lin         = this.Lin        (1:this.NAcq);
            this.Par         = this.Par        (1:this.NAcq);
            this.Sli         = this.Sli        (1:this.NAcq);
            this.Ave         = this.Ave        (1:this.NAcq);
            this.Phs         = this.Phs        (1:this.NAcq);
            this.Eco         = this.Eco        (1:this.NAcq);
            this.Rep         = this.Rep        (1:this.NAcq);
            this.Set         = this.Set        (1:this.NAcq);
            this.Seg         = this.Seg        (1:this.NAcq);
            this.Ida         = this.Ida        (1:this.NAcq);
            this.Idb         = this.Idb        (1:this.NAcq);
            this.Idc         = this.Idc        (1:this.NAcq);
            this.Idd         = this.Idd        (1:this.NAcq);
            this.Ide         = this.Ide        (1:this.NAcq);
            this.centerCol   = this.centerCol  (1:this.NAcq);
            this.centerLin   = this.centerLin  (1:this.NAcq);
            this.centerPar   = this.centerPar  (1:this.NAcq);
            this.IsReflected = this.IsReflected(1:this.NAcq);
            this.IsRawDataCorrect = this.IsRawDataCorrect(1:this.NAcq); %SRY;
            this.slicePos    = this.slicePos (:,1:this.NAcq);
            this.iceParam    = this.iceParam (:,1:this.NAcq);
            this.freeParam   = this.freeParam(:,1:this.NAcq);
            this.memPos      = this.memPos     (1:this.NAcq);
        
            this.NLin = max(this.Lin);
            this.NPar = max(this.Par);
            this.NSli = max(this.Sli);
            this.NAve = max(this.Ave);
            this.NPhs = max(this.Phs);
            this.NEco = max(this.Eco);
            this.NRep = max(this.Rep);
            this.NSet = max(this.Set);
            this.NSeg = max(this.Seg);
            this.NIda = max(this.Ida);
            this.NIdb = max(this.Idb);
            this.NIdc = max(this.Idc);
            this.NIdd = max(this.Idd);
            this.NIde = max(this.Ide);
            
            % ok, let us assume for now that all NCol and NCha entries are
            % the same for all mdhs:
            this.NCol = this.NCol(1);
            this.NCha = this.NCha(1);
                            
            this.dataDims = {'Col','Cha','Lin','Par','Sli','Ave','Phs',...
                'Eco','Rep','Set','Seg','Ida','Idb','Idc','Idd','Ide'};
            
            % to reduce the matrix sizes of non-image scans, the size 
            % of the refscan_obj()-matrix is reduced to the area of the
            % actually scanned acs lines (the outer part of k-space
            % that is not scanned is not filled with zeros)
            % this behaviour is controlled by flagSkipToFirstLine which is
            % set to true by default for everything but image scans
            if ~this.flagSkipToFirstLine
                % the output matrix should include all leading zeros
                this.skipLin = 0;
                this.skipPar = 0;
            else
                % otherwise, cut the matrix size to the start of the
                % first actually scanned line/partition (e.g. the acs/
                % phasecor data is only acquired in the k-space center)
                this.skipLin = min(this.Lin)-1;
                this.skipPar = min(this.Par)-1;
            end
            NLinAlloc = max(1, this.NLin - this.skipLin);
            NParAlloc = max(1, this.NPar - this.skipPar);

            this.fullSize = [ this.NCol this.NCha NLinAlloc NParAlloc...
                              this.NSli this.NAve this.NPhs this.NEco...
                              this.NRep this.NSet this.NSeg this.NIda...
                              this.NIdb this.NIdc this.NIdd this.NIde ];
                          
            this.dataSize = this.fullSize;

            if this.arg.removeOS
                this.dataSize(1) = this.NCol/2;
            end

            if this.arg.doAverage
                this.dataSize(6) = 1;
            end

            if this.arg.averageReps
                this.dataSize(9) = 1;
            end
            
            if this.arg.averageSets
                this.dataSize(10) = 1;
            end
        
            if this.arg.ignoreSeg
                this.dataSize(11) = 1;
            end

            % calculate sqzSize
            this.calcSqzSize;
            
            % calculate indices to target & source(raw)
            this.calcIndices;
            
            nByte = this.NCha*(this.freadInfo.szChannelHeader+8*this.NCol);
    
            % size for fread
            this.freadInfo.sz    = [2 nByte/8];
            % reshape size
            this.freadInfo.shape = [this.NCol+this.freadInfo.szChannelHeader/8 ...
                                   , this.NCha];
            % we need to cut MDHs from fread data
            this.freadInfo.cut   = this.freadInfo.szChannelHeader/8+1 ...
                              : this.NCol+this.freadInfo.szChannelHeader/8;
                          
                          
            % SRY: check that the number of raw data correction factors matches the
            % channel count              
            if (strcmp(this.softwareVersion, 'vb')) % not implemented/tested for vd, yet
                if (length(this.arg.rawDataCorrectionFactors) ~= this.NCha)
                    error('Number of raw data correction factors (%d) does not equal number of channels (%d)',...
                        length(rawDataCorrectionFactors), this.NCha);
                end
            end
           
        end
        
        
        function varargout = subsref(this, S)
            % this is where the magic happens
            % Now seriously. Overloading of the subsref-method and working 
            % with a gazillion indices got really messy really fast. At
            % some point, I should probably clean this code up a bit. But
            % good news everyone: It seems to work.
            
            switch S(1).type
                case '.' 
                    % We don't want to manage method/variable calls, so we'll 
                    % simply call the built-in subsref-function in this case.
                    % Note, that this has the limitation that there's no way
                    % to find out whether the original call was terminated
                    % with a semicolon! So "this.that" won't produce any
                    % output and is identical to "this.that;"
                    if nargout == 0
                        builtin('subsref', this, S);
                    else
                        varargout      = cell(1, nargout);
                        [varargout{:}] = builtin('subsref', this, S);
                    end
                    return;
                case {'()','{}'}
                otherwise
                    error('operator not supported');
            end
           
            [selRange selRangeSz outSize] = this.calcRange(S(1));
            
            tmp = reshape(1:prod(double(this.fullSize(3:end))), this.fullSize(3:end));
            tmp = tmp(selRange{3:end});
            cIxToRaw = this.ixToRaw(tmp); clear tmp;
            cIxToRaw = cIxToRaw(:);
            % delete all entries that point to zero (the "NULL"-pointer)
            notAcquired = (cIxToRaw == 0);
            cIxToRaw (notAcquired) = []; clear notAcquired;
            
            % calculate cIxToTarg for possibly smaller, shifted + segmented
            % target matrix:
            cIx = zeros(14, numel(cIxToRaw), 'single');
            cIx( 1,:) = this.Lin(cIxToRaw) - this.skipLin;
            cIx( 2,:) = this.Par(cIxToRaw) - this.skipPar;
            cIx( 3,:) = this.Sli(cIxToRaw);
            if this.arg.doAverage
                cIx( 4,:) = 1;
            else
                cIx( 4,:) = this.Ave(cIxToRaw);
            end
            cIx( 5,:) = this.Phs(cIxToRaw);
            cIx( 6,:) = this.Eco(cIxToRaw);
            if this.arg.averageReps
                cIx( 7,:) = 1;
            else
                cIx( 7,:) = this.Rep(cIxToRaw);
            end
            if this.arg.averageSets
                cIx( 8,:) = 1;
            else
                cIx( 8,:) = this.Set(cIxToRaw);
            end
            if this.arg.ignoreSeg
                cIx( 9,:) = 1;
            else
                cIx( 9,:) = this.Seg(cIxToRaw);
            end
            cIx(10,:) = this.Ida(cIxToRaw);
            cIx(11,:) = this.Idb(cIxToRaw);
            cIx(12,:) = this.Idc(cIxToRaw);
            cIx(13,:) = this.Idd(cIxToRaw);
            cIx(14,:) = this.Ide(cIxToRaw);
            
            % make sure that indices fit inside selection range
            for k=3:numel(selRange)
                tmp = cIx(k-2,:);
                for l=1:numel(selRange{k})
                    cIx(k-2,tmp==selRange{k}(l)) = l;
                end
            end

            cIxToTarg = sub2ind_double(selRangeSz(3:end),cIx(1,:),cIx(2,:),cIx(3,:),...
                cIx(4,:),cIx(5,:),cIx(6,:),cIx(7,:),cIx(8,:),cIx(9,:),...
                cIx(10,:),cIx(11,:),cIx(12,:),cIx(13,:),cIx(14,:));
            
            mem = this.memPos(cIxToRaw);
            % sort mem for quicker access, sort cIxToTarg/Raw accordingly
            [mem,ix]  = sort(mem);
            cIxToTarg = cIxToTarg(ix);
            cIxToRaw  = cIxToRaw(ix);
            clear ix;

            out = complex(zeros(outSize,'single'));
            out = reshape(out, selRangeSz(1), selRangeSz(2), []);

            % counter for proper scaling of averages/segments
            count_ave = zeros([1 1 size(out,3)],'single');
            
            % subsref overloading makes this.that-calls slow, so we need to
            % avoid them whenever possible
            szScanHeader = this.freadInfo.szScanHeader;
            readSize     = this.freadInfo.sz;
            readShape    = this.freadInfo.shape;
            readCut      = this.freadInfo.cut;
            cutOS        = this.NCol/4+1:this.NCol*3/4;
            bRemoveOS    = this.arg.removeOS;
            bIsReflected = this.IsReflected;
            %SRY store information about raw data correction
            bDoRawDataCorrect = this.arg.doRawDataCorrect;
            bIsRawDataCorrect = this.IsRawDataCorrect;
            if (bDoRawDataCorrect)
                rawDataCorrect = this.arg.rawDataCorrectionFactors;
            end
            
            fid = fopen(this.filename);
            
            for k=1:numel(mem)
                % skip scan header            
                fseek(fid,mem(k) + szScanHeader,'bof');
                raw = fread(fid, readSize, 'float=>single').';
                raw = reshape( complex(raw(:,1), raw(:,2)), readShape);
                raw = raw(readCut,:);

                %SRY apply raw data correction if necessary
                if ( bDoRawDataCorrect && bIsRawDataCorrect(cIxToRaw(k)) )
                    %there are two ways to do this: multiply where the flag is
                    %set, or divide where it is not set.  There are significantly
                    %more points without the flag, so multiplying is more
                    %efficient
                    raw = bsxfun(@times, raw, rawDataCorrect);
                end
            
                % select channels
                raw = raw(:,selRange{2});
                
                if bRemoveOS
                    % remove oversampling in read
                    raw          = ifft(raw);
                    raw(cutOS,:) = [];
                    raw          = fft(raw);
                end

                if bIsReflected(cIxToRaw(k))
                    raw = raw(end:-1:1,:);
                end
                
                % select columns and sort data
                out(:,:,cIxToTarg(k))       = out(:,:,cIxToTarg(k)) +...
                                              raw(selRange{1},:);

                count_ave(1,1,cIxToTarg(k)) = count_ave(1,1,cIxToTarg(k)) + 1;
            end    
            fclose(fid);
                    
            % proper scaling (we don't want to sum our data but average it)
            if ~this.flagNoWeightedAverage
                count_ave = 1./max(1,count_ave);
                out       = bsxfun(@times,out,count_ave);
            end
            
            out = reshape(out,outSize);
 
            % For a call of type data{:,:,1:3} matlab expects more than one
            % output variable (three in this case) and will throw an error 
            % otherwise. This is a lazy way (and the only one I know of) to 
            % fix this.
            varargout    = cell(1, nargout);
            varargout{1} = out;
        end
        

        function unsorted = unsorted(this)
            % returns the unsorted data [NCol,NCha,#samples in acq. order]
            szScanHeader = this.freadInfo.szScanHeader;
            readSize     = this.freadInfo.sz;
            readShape    = this.freadInfo.shape;
            readCut      = this.freadInfo.cut;
            cutOS        = this.NCol/4+1:this.NCol*3/4;
            bRemoveOS    = this.arg.removeOS;
            bIsReflected = this.IsReflected;
            %SRY store information about raw data correction
            bDoRawDataCorrect = this.arg.doRawDataCorrect;
            bIsRawDataCorrect = this.IsRawDataCorrect;
            if (bDoRawDataCorrect)
                rawDataCorrect = this.arg.rawDataCorrectionFactors;
            end
            mem = this.memPos;
            unsorted = complex(zeros([this.dataSize(1) this.NCha this.NAcq],'single'));
            fid = fopen(this.filename);
            for k=1:this.NAcq
                % skip scan header            
                fseek(fid,mem(k) + szScanHeader,'bof');
                raw = fread(fid, readSize, 'float=>single').';
                raw = reshape( complex(raw(:,1), raw(:,2)), readShape);
                raw = raw(readCut,:);
                if ( bDoRawDataCorrect && bIsRawDataCorrect(k) )
                    %SRY apply raw data correction if necessary
                    raw = bsxfun(@times, raw, rawDataCorrect);
                end
                if bRemoveOS
                    % remove oversampling in read
                    raw          = ifft(raw);
                    raw(cutOS,:) = [];
                    raw          = fft(raw);
                end
                if bIsReflected(k)
                    raw = raw(end:-1:1,:);
                end
                unsorted(:,:,k) = raw;
            end
            fclose(fid);
        end % of get.unsorted()
        
        
        function set.flagRemoveOS(this,val)
            % set method for removeOS 
            this.arg.removeOS = logical(val);
            
            % we also need to recalculate our data size:
            if this.arg.removeOS
                this.dataSize(1) = this.NCol(1)/2;
                this.sqzSize(1)  = this.NCol(1)/2;
            else
                this.dataSize(1) = this.NCol(1);
                this.sqzSize(1)  = this.NCol(1);
            end
        end
        
        
        function out = get.flagRemoveOS(this)
            out = this.arg.removeOS;
        end
        
        
        function set.flagDoAverage(this,val)
            % set method for doAverage 
            this.arg.doAverage = logical(val);
            
            if this.arg.doAverage
                this.dataSize(6) = 1;
            else
                this.dataSize(6) = this.NAve;
            end
            
            % update sqzSize
            this.calcSqzSize;
        end
        
        
        function out = get.flagDoAverage(this)
            out = this.arg.doAverage;
        end
        
        function set.flagAverageReps(this,val)
            % set method for doAverage 
            this.arg.averageReps = logical(val);
            
            if this.arg.averageReps
                this.dataSize(9) = 1;
            else
                this.dataSize(9) = this.NRep;
            end
            
            % update sqzSize
            this.calcSqzSize;
        end
        
        
        function out = get.flagAverageReps(this)
            out = this.arg.averageReps;
        end
        
        
        function set.flagAverageSets(this,val)
            % set method for doAverage 
            this.arg.averageSets = logical(val);
            
            if this.arg.averageSets
                this.dataSize(10) = 1;
            else
                this.dataSize(10) = this.NSet;
            end
            
            % update sqzSize
            this.calcSqzSize;
        end
        
        
        function out = get.flagAverageSets(this)
            out = this.arg.averageSets;
        end
        
        
        function set.flagSkipToFirstLine(this,val)
            val = logical(val);
            if val ~= this.arg.skipToFirstLine
                this.arg.skipToFirstLine = val;

                if this.arg.skipToFirstLine
                    this.skipLin = min(this.Lin)-1;
                    this.skipPar = min(this.Par)-1;
                else
                    this.skipLin = 0;
                    this.skipPar = 0;
                end
                NLinAlloc = max(1, this.NLin - this.skipLin);
                NParAlloc = max(1, this.NPar - this.skipPar);
                this.fullSize(3:4) = [NLinAlloc NParAlloc];
                this.dataSize(3:4) = this.fullSize(3:4);

                % update sqzSize
                this.calcSqzSize;
                % update indices
                this.calcIndices;
            end
        end
        
        
        function out = get.flagSkipToFirstLine(this)
            out = this.arg.skipToFirstLine;
        end
        
        
        function set.flagIgnoreSeg(this,val)
            % set method for ignoreSeg 
            this.arg.ignoreSeg = logical(val);
            
            if this.arg.ignoreSeg
                this.dataSize(11) = 1;
            else
                this.dataSize(11) = this.NSeg;
            end
            
            % update sqzSize
            this.calcSqzSize;
        end
        
        
        function out = get.flagIgnoreSeg(this)
            out = this.arg.ignoreSeg;
        end
        
        %SRY: accessor methods for raw data correction
        function out = get.flagDoRawDataCorrect(this)
            out = this.arg.doRawDataCorrect;
        end
      
        function set.flagDoRawDataCorrect(this, val)
            val = logical(val);
            if (val == true && strcmp(this.softwareVersion, 'vd'))
                error('raw data correction for VD not supported/tested yet');
            end
               
            this.arg.doRawDataCorrect = val;
        end
      
        function out = get.RawDataCorrectionFactors(this)
            out = this.arg.rawDataCorrectionFactors;
        end
        
        function set.RawDataCorrectionFactors(this, val)
            %this may not work if trying to set the factors before NCha has
            %a meaningful value (ie before calling clean)
            if (~isrow(val) || length(val) ~= this.NCha)
                error('RawDataCorrectionFactors must be a 1xNCha row vector');
            end
            this.arg.rawDataCorrectionFactors = val;
        end
      
    end
    
    methods (Access='protected')
        % helper functions
        
        function [selRange selRangeSz outSize] = calcRange(this,S)
            
            switch S.type
                case '()'
                    bSqueeze = false;
                case '{}'
                    bSqueeze = true;
            end
                    
            selRange = num2cell(ones(1,numel(this.dataSize)));
            outSize  = ones(1,numel(this.dataSize));
            
            if ( isempty(S.subs) || strcmpi(S.subs(1),'') )
                % obj(): shortcut to select all data
                % unfortunately, matlab does not allow the statement
                % obj{}, so we can't use it...
                % alternative: obj{''} (obj('') also works)
                for k=1:numel(this.dataSize)
                    selRange{k}   = 1:this.dataSize(k);
                end
                if ~bSqueeze
                    outSize = this.dataSize;
                else
                    outSize = this.sqzSize;
                end
            else    
                for k=1:numel(S.subs)
                    if ~bSqueeze
                        cDim = k; % nothing to do
                    else
                        % we need to rearrange selRange from squeezed
                        % to original order
                        cDim = find(strcmp(this.dataDims,this.sqzDims{k}) == 1);
                    end
                    if strcmp(S.subs{k},':')
                        if k<numel(S.subs)
                            selRange  {cDim} = 1:this.dataSize(cDim);
                        else % all later dimensions selected and 'vectorized'!
                            for l=cDim:numel(this.dataSize)
                                selRange{l} = 1:this.dataSize(l);
                            end
                            outSize(k) = prod(double(this.dataSize(cDim:end)));
                            break; % jump out ouf for-loop
                        end
                    elseif isnumeric(S.subs{k})
                        selRange{cDim} = single(S.subs{k});
                    else
                        error('unknown string in brackets (e.g. 1:end does not work here)');
                    end
                    outSize(k) = numel(selRange{cDim});
                end
            end
            
            for k=1:numel(selRange)
               if max(selRange{k}) > this.dataSize(k)
                   error('selection out of range');
               end
            end
            
            selRangeSz = ones(1,numel(this.dataSize));
            for k=1:numel(selRange)
                selRangeSz(k) = numel(selRange{k});
            end          
            
            % now select all averages in case doAverage is selected
            if this.arg.doAverage
                selRange{6}  = 1:this.fullSize(6);
            end
            % now select all repetitions in case averageReps is selected
            if this.arg.averageReps
                selRange{9}  = 1:this.fullSize(9);
            end
            % now select all sets in case averageSets is selected
            if this.arg.averageSets
                selRange{10}  = 1:this.fullSize(10);
            end
            % now select all segments in case ignoreSeg is selected
            if this.arg.ignoreSeg
                selRange{11} = 1:this.fullSize(11);
            end
            
        end
        
        function calcSqzSize(this)
            % calculate sqzSize and sqzDims
            this.sqzSize    = [];
            this.sqzDims    = [];
            this.sqzSize(1) = this.dataSize(1);
            this.sqzDims{1} = 'Col';
            c = 1;
            for k=2:numel(this.dataSize)
                if this.dataSize(k)>1
                    c = c+1;
                    this.sqzSize(c) = this.dataSize(k);
                    this.sqzDims{c} = this.dataDims{k};
                end
            end
        end
        
        function calcIndices(this)
            % calculate indices to target & source(raw)
            LinIx     = this.Lin  - this.skipLin;
            ParIx     = this.Par  - this.skipPar;
            this.ixToTarget = sub2ind_double(this.fullSize(3:end),...
                LinIx, ParIx, this.Sli, this.Ave, this.Phs, this.Eco,...
                this.Rep, this.Set, this.Seg, this.Ida, this.Idb,...
                this.Idc, this.Idd, this.Ide);
            
            % now calculate inverse index
            % inverse index of lines that are not measured is zero
            this.ixToRaw = zeros(1,prod(this.fullSize(3:end)),'double');

            % subsref overloading makes this.that-calls slow, so we need to
            % avoid them whenever possible
            ixToTarg = this.ixToTarget;
            for k=1:numel(ixToTarg)
                this.ixToRaw(ixToTarg(k)) = k;
            end
        end

    end
    
end


%%%%%%%%%%%% helper functions %%%%%%%%%%%
function ndx = sub2ind_double(sz,varargin)
%SUB2IND_double Linear index from multiple subscripts.
%   Works like sub2ind but always returns double
%   also slightly faster, but no checks
%========================================
    sz  = double(sz);
    ndx = double(varargin{end}) - 1;
    for i = length(sz)-1:-1:1
        ix  = double(varargin{i});
        ndx = sz(i)*ndx + ix-1;
    end
    ndx = ndx + 1;
end
