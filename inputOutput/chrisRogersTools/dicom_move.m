% Advance to the point in the file where the target group and
% element are and return the length of the field.
% Orginally:
% Greg Reynolds 26-January-2005.
% Vastly updated:
% Greg Reynolds 17-June-2005.
%
% Now with in-sequence data support.

function length = dicom_move(fd, strGroup, strElement)

fprintf('\nSearching for target element (%s, %s)...', strGroup, strElement');

% these are all the VRs that have length 2
VR_short_length = struct('strings', ...
    {'AE', 'AS', 'AT', 'CS', 'DA', 'DS', 'DT', ...
        'FL', 'FD', 'IS', 'LO', 'LT', 'OF', 'PN', 'SH', ...
        'SL', 'ST', 'SS', 'TM', 'UI', 'UL', 'US' });

dims = size(VR_short_length);
number_of_short_vrs = dims(2);

% these are all the VRs so that we can establish implicit VRs
% without lots of DICOM knowledge
VRs = struct('strings', ...
    {'AE', 'AS', 'AT','CS', 'DA', 'DS', 'DT', 'FL', ...
        'FD', 'IS', 'LO', 'LT', 'OB', 'OF', 'OW', 'PN', ...
        'SH', 'SL', 'SQ', 'ST', 'SS', 'TM', 'UI', 'UL', ...
        'UN', 'US', 'UT'});

dims = size(VRs);
number_of_vrs = dims(2);

done = 0;

while ~done,   
   
   current_tag = fread(fd, 2, 'uint16','l');
   current_vr = fread(fd, 2, 'schar','l');   
    
   if feof(fd)
       done = 1;
       length = 0;
       fprintf('\nReached end of file without match.');
       break;
   else               
        strGroupCurrent = sprintf('%X', current_tag(1));
        strElementCurrent = sprintf('%X', current_tag(2));
        strVRCurrent = sprintf('%c', current_vr);

        % first of all, check with this is an implicit VR
        explicit_vr = 0;
        for n = 1:1:number_of_vrs
            if strcmp(VRs(n).strings, strVRCurrent)
                explicit_vr = 1;
                break;
            end
        end
        
        % it was an implicit VR
        if explicit_vr == 0
            
            % adjust the file pointer back the two-bytes we tentatively
            % read in as being the VR
            fseek(fd, -2, 'cof');
            
            % possibly need to read in zero padding here?
            current_length = fread(fd, 1, 'uint32','l');

            % if the length is undefined, just drop out and move
            % to next element...
            if ~strcmp(sprintf('%X', current_length), 'FFFFFFFF')

                if strcmp(strGroupCurrent, strGroup)
                    if strcmp(strElementCurrent, strElement)
                        length = current_length;
                        done = 1;
                        break;
                    end
                end
                
                if done == 0
                    fread(fd, current_length, 'uchar','l');      
                end
            end
            
        % it was an explicit VR
        else                       
            size_length = 4;        
            
            % check to see whether it has a short length
            for n = 1:1:number_of_short_vrs
                if strcmp(VR_short_length(n).strings, strVRCurrent)
                    size_length = 2;
                end
            end
            
            % note that implicit VRs always have 32-bit length
            if size_length == 2
                current_length = fread(fd, 1, 'uint16','l');
            else
                zeropadding = fread(fd, 2, 'uchar','l');
                current_length = fread(fd, 1, 'uint32','l');	    
            end           
            
            % now see if this was the field we wanted
            if strcmp(strGroupCurrent, strGroup)
                if strcmp(strElementCurrent, strElement)
                    length = current_length;
                    done = 1;
                    break;
                end
            end
                                   
            % some implicit VRs, e.g. an SQ can have undefined
            % length if they have 32-bit length, so if that
            % isn't the case, proceed and read, otherwise just
            % skip to the next element            
            if ~strcmp(sprintf('%X', current_length), 'FFFFFFFF')   
                % it wasn't, so advance pointer
                if done == 0
                    fread(fd, current_length, 'uchar','l');
                end
            end	
        end    
    end
end


