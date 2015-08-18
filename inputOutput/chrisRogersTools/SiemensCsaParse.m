function [dicomInfo] = SiemensCsaParse(dicomInfo)
% Read values from the 'SIEMENS CSA HEADER' private group in a
% DICOM file produced by a Siemens MR scanner.
%
% [dicomInfo] = SiemensCsaParse(dicomInfo)
%
% dicomInfo is a Siemens MR DICOM file loaded using Matlab's dicominfo()
%        OR a string specifying the filename of a DICOM file to load.
%
% output is the original structure with an extra field 'csa' that contains
% the decoded 'SIEMENS CSA HEADER' block

% Inspiration for this comes from the gdcm tool (http://gdcm.sf.net/),
% specifically gdcmCSAHeader.cxx lines 490-65.

% Possible TODO:
% 1. Migrate to C++ MEX file.
% 2. Investigate whether there's a Siemens product .dll file (e.g.
%    Z:\n4\pkg\MrServers\MrProtSrv) for parsing these headers.
% 3. Understand meaning of syngodt (Syngo Datatype ???) field.

% Copyright: Chris Rodgers (University of Oxford), 2008-11.
% All rights reserved.

% $Id: SiemensCsaParse.m 5022 2012-02-01 11:26:39Z will $

if ischar(dicomInfo)
    dicomInfo = dicominfo(dicomInfo);
end

%% Locate SIEMENS CSA HEADER
tagId = getDicomPrivateTag(dicomInfo,'0029','SIEMENS CSA HEADER');

if isempty(tagId)
    % Give warning and return with unmodified dicomInfo
    warning('RodgersSpectroTools:NoSiemensCsaHeader','This DICOM file does not contain a ''SIEMENS CSA HEADER''.')
    return
end

%% Process the two parameter fields
SiemensCsaParse_ReadDicomTag(['Private_0029_' tagId '10']);
SiemensCsaParse_ReadDicomTag(['Private_0029_' tagId '20']);

return





%% Main code
function SiemensCsaParse_ReadDicomTag(strTag)
    currdx=0;
    
    if ~strcmp(char(private_read(4)),'SV10') || ~all(private_read(4)==[4 3 2 1])
        error('Unsupported CSA block format');
    end
    
    % This parsing code is translated from gdcm (http://gdcm.sf.net/)
    numElements = double(private_readuint32(1));
    
    % Sanity check
    if private_readuint32(1)~=77
        error('Unsupported CSA block format');
    end
    
    for tagdx=1:numElements
        tagName = private_readstring(64);
        
        % Fix up tagName
        tagName(tagName == '-') = [];
        
        vm = private_readuint32(1);
        vr = private_readstring(4);
        syngodt = private_readuint32(1);
        nitems = double(private_readuint32(1));
        
        checkbit = private_readuint32(1);
        
        if checkbit ~= 77 && checkbit ~= 205
            error('Unsupported CSA block format');
        end
        
        data = {};
        for itemdx=1:nitems
            header = double(private_readuint32(4));
            
            if (header(3) ~= 77 && header(3) ~= 205) || ...
                    (header(1) ~= header(2)) || ...
                    (header(1) ~= header(4))
                error('Unsupported CSA block format');
            end
            
            data{itemdx} = private_readstring(header(1));
            
            % Dump junk up to DWORD boundary
            private_read(mod(mod(4-header(1),4),4));
        end
        
        % Store this in the csa structure
        switch vr
            case {'CS', 'LO', 'LT', 'SH', 'SS', 'UI', 'UT', 'UN'} % Strings and unknown byte string
                if numel(data) < vm
                    % Pad if necessary. Siemens CSA format omits null strings.
                    data{vm} = '';
                end
                
                if vm == 1
                    dicomInfo.csa.(tagName) = data{1};
                else
                    dicomInfo.csa.(tagName) = data(1:vm);
                end
            case {'DS', 'FD', 'FL', 'IS', 'SL', 'ST', 'UL', 'US'} % Numbers
                dataNumeric = arrayfun(@str2double,data);
                
                if numel(dataNumeric) < vm
                    % Zero pad if necessary. Siemens CSA format omits zeros.
                    dataNumeric(vm) = 0;
                end
                
                dicomInfo.csa.(tagName) = dataNumeric(1:vm);
            otherwise
                warning('RodgersSpectroTools:UnknownVrType','Unknown VR type: "%s".',vr)
        end
    end
        
        
    %% Helper functions to simulate file I/O
    function [out] = private_read(numBytes)
        out = dicomInfo.(strTag)(currdx+(1:numBytes)).';
        currdx=currdx+numBytes;
    end
    
    function [out] = private_readuint32(num)
        out=typecast(private_read(4*num),'uint32');
    end
    
    function [out] = private_readstring(maxchar)
        out = reshape(char(private_read(maxchar)),1,[]);
        terminator = find(out==0,1);
        if numel(terminator)>0
            out=out(1:(terminator-1));
        end
    end
    
end

end
