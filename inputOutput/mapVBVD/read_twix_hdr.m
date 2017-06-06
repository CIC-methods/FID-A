function [prot,rstraj] = read_twix_hdr(fid)

% function to read raw data header information from siemens MRI scanners 
% (currently VB and VD software versions are supported and tested).
%
% Author: Philipp Ehses (philipp.ehses@tuebingen.mpg.de), Mar/11/2014
      
    nbuffers = fread(fid, 1,'uint32');
    
    prot = [];
    for b=1:nbuffers
        %now read string up to null termination     
        bufname = fread(fid, 10, 'uint8=>char').';
        bufname = regexp(bufname, '^\w*', 'match');
        bufname = bufname{1};
        fseek(fid, numel(bufname)-9, 'cof');        
        buflen         = fread(fid, 1,'uint32');
        buffer         = fread(fid, buflen, 'uint8=>char').';
        buffer         = regexprep(buffer,'\n\s*\n',''); % delete empty lines
        prot.(bufname) = parse_buffer(buffer);
    end
    
    if nargout>1
        rstraj = [];
        if isfield(prot.Meas,'alRegridMode') && prot.Meas.alRegridMode(1)>1
            ncol           = prot.Meas.alRegridDestSamples(1);
            dwelltime      = prot.Meas.aflRegridADCDuration(1)/ncol;
            gr_adc         = zeros(1,ncol,'single');
            start          = prot.Meas.alRegridRampupTime(1) - (prot.Meas.aflRegridADCDuration(1)-prot.Meas.alRegridFlattopTime(1))/2;
            time_adc       = start + dwelltime * (0.5:ncol);            
            ixUp           = time_adc <= prot.Meas.alRegridRampupTime(1);
            ixFlat         = (time_adc <= prot.Meas.alRegridRampupTime(1)+prot.Meas.alRegridFlattopTime(1)) & ~ixUp;
            ixDn           = ~ixUp & ~ixFlat;
            gr_adc(ixFlat) = 1;
            if prot.Meas.alRegridMode(1) == 2  
                % trapezoidal gradient
                gr_adc(ixUp)   = time_adc(ixUp)/prot.Meas.alRegridRampupTime(1);
                gr_adc(ixDn)   = 1 - (time_adc(ixDn)-prot.Meas.alRegridRampupTime(1)-prot.Meas.alRegridFlattopTime(1))/prot.Meas.alRegridRampdownTime(1);
            elseif prot.Meas.alRegridMode(1) == 4  
                % sinusoidal gradient
                gr_adc(ixUp)   = sin(pi/2*time_adc(ixUp)/prot.Meas.alRegridRampupTime(1));
                gr_adc(ixDn)   = sin(pi/2*(1+(time_adc(ixDn)-prot.Meas.alRegridRampupTime(1)-prot.Meas.alRegridFlattopTime(1))/prot.Meas.alRegridRampdownTime(1)));
            else
                warning('regridding mode unknown');
                return;
            end
            % make sure that gr_adc is always positive (rstraj needs to be
            % strictly monotonic):
            gr_adc = max(gr_adc,1e-4);
            rstraj = (cumsum(gr_adc(:)) - ncol/2)/sum(gr_adc(:));
            rstraj = rstraj - rstraj(ncol/2+1);
            % scale rstraj by kmax (only works if all slices have same FoV!!!)
            kmax = prot.MeasYaps.sKSpace.lBaseResolution/...
                prot.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV;
            rstraj = kmax * rstraj;
        end
    end
    
end        


function prot = parse_buffer(buffer)
    [ascconv, xprot] = regexp(buffer,'### ASCCONV BEGIN[^\n]*\n(.*)\s### ASCCONV END ###','tokens','split');

    if ~isempty(ascconv)
        ascconv = ascconv{:}{:};
        prot = parse_ascconv(ascconv);
    else
        prot = struct();
    end

    if ~isempty(xprot)
        xprot = strcat(xprot{:}); % bug fix by Qiuting Wen (added strcat)
        xprot = parse_xprot(xprot);
        if isstruct(xprot)
            name   = cat(1,fieldnames(prot),fieldnames(xprot));
            val    = cat(1,struct2cell(prot),struct2cell(xprot));
            [~,ix] = unique(name);
            prot   = cell2struct(val(ix),name(ix));
        end
    end
end


function xprot = parse_xprot(buffer)
    xprot = [];
    tokens = regexp(buffer, '<Param(?:Bool|Long|String)\."(\w+)">\s*{([^}]*)','tokens');
    tokens = [tokens, regexp(buffer, '<ParamDouble\."(\w+)">\s*{\s*(<Precision>\s*[0-9]*)?\s*([^}]*)','tokens')];
    for m=1:numel(tokens)
        name         = char(tokens{m}(1));
        % field name has to start with letter
        if (~isletter(name(1)))
            name = strcat('x', name);
        end

        value = char(strtrim(regexprep(tokens{m}(end), '("*)|( *<\w*> *[^\n]*)', '')));
        value = regexprep(value, '\s*', ' ');

        try %#ok<TRYNC>
            value = eval(['[' value ']']);  % inlined str2num()
        end

        xprot.(name) = value;
    end
end


function mrprot = parse_ascconv(buffer)  
    mrprot = [];    
    % [mv] was: vararray = regexp(buffer,'(?<name>\S*)\s*=\s(?<value>\S*)','names');
    vararray = regexp(buffer,'(?<name>\S*)\s*=\s*(?<value>\S*)','names');
    
    isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    
    for var=vararray

        try
            value = eval(['[' var.value ']']);  % inlined str2num()
        catch
            value = var.value;
        end
        
        % now split array name and index (if present)
        v = regexp(var.name,'(?<name>\w*)\[(?<ix>[0-9]*)\]|(?<name>\w*)','names');

        cnt = 0;
        tmp = cell(2, numel(v));

        breaked = false;
        for k=1:numel(v)
            if isOctave
                vk = v{k};
                if iscell(vk.name)
                    % lazy fix that throws some info away
                    vk.name = vk.name{1};
                    vk.ix   = vk.ix{1};
                end
            else
                vk = v(k);
            end
            if ~isletter(vk.name(1))
                breaked = true;
                break;
            end
            cnt = cnt+1;
            tmp{1,cnt} = '.';
            tmp{2,cnt} = vk.name;

            if ~isempty(vk.ix)
                cnt = cnt+1;
                tmp{1,cnt} = '{}';
                tmp{2,cnt}{1} = 1 + str2double(vk.ix);
            end
        end
        if ~breaked && ~isempty(tmp)
            S = substruct(tmp{:});
            mrprot = subsasgn(mrprot,S,value);
        end
    end 
end
