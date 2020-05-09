function ParList = read_IceProgramParam(file_path)
%
% read_ascconv Read ascconv header part of DICOM and Siemens raw data
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function cuts out the ascconv header part of DICOM and Siemens raw data and searches for Parameters within this header. These
%
%
% [ParList,ascconv] = read_ascconv(file_path)
%
% Input: 
% -         file_path                     ...     Path of file.
%
% Output:
% -         ParList                       ...     Structure giving all the Parameters. It contains among many others:
%
%           -- ParList.total_channel_no_measured         - Number of receive-channels that were measured            
%           -- ParList.total_channel_no_reco             - Number of receive-channels that were are in the file (DICOM)             
%           -- ParList.Dwelltimes                        - Dwelltime             
%           -- ParList.LarmorFreq                        - LarmorFrequency                
%           -- ParList.ThreeD_flag                       - flag telling you if measurement is real 3D measurement or 2D/Multislice                     
%           -- ParList.AsymmetricEcho                    - If AssymetricEcho was allowed                    
%           -- ParList.InterleavedSliceAcquisition       - If Slices were not measured consecutively like 1,2,3,4,... but e.g. 1,3,2,4         
%           -- ParList.nFreqEnc                          - Number of measured points in frequency encoding direction       
%           -- ParList.nPhasEnc                          - Number of measured points in phase encoding direction           
%           -- ParList.nPartEnc                          - Number of measured points in partition encoding direction (only for 3D)
%           ...
%           ...
%           ...
%
%
% -         ascconv                       ...     cell array of the ascconv header: ascconv{:,1} = ParameterName, ascconv{:,2} = ParameterValue 
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None






%% 0. Preparations


% Define for which entries the ascconv should be searched for
% Search for these entries in the ascconv header part:


% open file
fid = fopen(file_path,'r');
headersize = fread(fid,1, 'uint32');


%% 1. Track down & save element: ASCCONV

begin_found = 0;
sLine = 0;
ParList = [];
CurFind = 0;
while(sLine > -1)
    
    sLine = fgets(fid); % get string line
    CurFilePos_Bytes = ftell(fid);
    if(CurFilePos_Bytes > headersize)
        break
    end
    
    
        
        
        if(not(~contains(sLine,'alICEProgramPara')))
            CurFind = CurFind + 1;
            ParList{CurFind} = sLine;
            begin_found = true;                                    
        else
            continue                                                % If current line is not begin of ascconv --> read next line
        end
        
     
        
end



%% 2. Display error & stop if no Ascconv found

if(not(begin_found))
    display(['Pfui Toifel! You gave me a file without any ascconv, I cannot digest that! Please remember that I am NOT an omnivore.' char(10) ...
             'I will stop here . . .'])
    return
end


%% 3. Convert ParList

% For now only extract the values, not their names...
ParList = regexp(ParList{1}, '{.*}','match');
ParList = regexprep(ParList{1}, '{','[');
ParList = regexprep(ParList, '}',']');
ParList2.Values = eval(ParList);
ParList = ParList2; clear ParList2





%% 5. Postparations

fclose(fid);

end





