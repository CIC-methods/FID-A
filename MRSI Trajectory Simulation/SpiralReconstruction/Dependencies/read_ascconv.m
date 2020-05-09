function [ParList,ascconv] = read_ascconv(file_path)
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
ParList_Search =  { ...
'asCoilSelectMeas\[0\]\.asList\[\d+\]\.lRxChannelConnected',	...     % 1
'Coil_Dumms',                                                   ...     % 27    Just a Dummy to have the parameter in the right place  
'sRXSPEC\.alDwellTime\[\d+\]',                                  ...     % 18
'sTXSPEC\.asNucleusInfo\[0\]\.lFrequency',                      ...     % 19	Only take the first one [0], why should there be more?
'sKSpace\.ucDimension',                                         ...     % 8     Is this the Parameter that is different for 3D vs. 2D acquisitions????? 
'sKSpace\.ucAsymmetricEchoAllowed',                             ...     % 11
'sKSpace\.ucMultiSliceMode',                                    ...     % 12   
'sKSpace\.lBaseResolution',                                     ...     % 2
'sKSpace\.lPhaseEncodingLines',                                 ...     % 3
'sKSpace\.lPartitions',                                         ...     % 6  
'sSpecPara\.lFinalMatrixSizeRead',                              ...     % 4  
'sSpecPara\.lFinalMatrixSizePhase',                             ...     % 5
'sSpecPara\.lFinalMatrixSizeSlice',                             ...     % 34    
'sSpecPara\.lVectorSize',                                       ...     % 9
'sSpecPara\.ucRemoveOversampling',                              ...     % 10
'sSpecPara.sVoI\.dPhaseFOV',                                    ...     % 23
'sSpecPara.sVoI\.dReadoutFOV',                                  ...     % 24
'sSpecPara.sVoI\.dThickness',                                   ...     % 25    
'sSliceArray\.lSize',                                           ...     % 7  
'sSliceArray.asSlice\[\d+\].dThickness',                        ...     % 20 
'sSliceArray\.asSlice\[\d+\]\.dPhaseFOV',                       ...     % 13
'sSliceArray\.asSlice\[\d+\]\.dReadoutFOV',                     ...     % 14
'FOV_Partition_dumms',                                          ...     % 22	Just a Dummy to have the parameter in the right place
'sSpecPara\.sVoI\.sPosition\.dSag',								...     % 43    Sagittal = forehead-backhead direction = "Read"-direction (normally)
'sSpecPara\.sVoI\.sPosition\.dCor',								...     % 44    Coronal = left-right directions = "Phase"-direction (normally)
'sSpecPara\.sVoI\.sPosition\.dTra',								...     % 45    Transversal = up-down direction = "Partition or Slice"-direction (normally)
'sSpecPara\.sVoI\.sNormal\.dSag',								...     % 46    The x-component of the normal vector of the measured slice
'sSpecPara\.sVoI\.sNormal\.dCor',								...     % 47    y-component
'sSpecPara\.sVoI\.sNormal\.dTra',								...     % 48    z-component
'sSpecPara\.sVoI\.dInPlaneRot',									...     % 49    The InPlane (InSlice) rotation, so the rotation around the normal vector given by the upper three components
'sSliceArray\.asSlice\[\d+\]\.sPosition\.dSag',                 ...     % 15    Sagittal = forehead-backhead direction = "Read"-direction (normally)
'sSliceArray\.asSlice\[\d+\]\.sPosition\.dCor',                 ...     % 16    Coronal = left-right directions = "Phase"-direction (normally)
'sSliceArray\.asSlice\[\d+\]\.sPosition\.dTra',                 ...     % 17    Transversal = up-down direction = "Partition or Slice"-direction (normally)
'sSliceArray\.asSlice\[\d+\]\.sNormal\.dSag',                   ...     % 28    The x-component of the normal vector of the measured slice
'sSliceArray\.asSlice\[\d+\]\.sNormal\.dCor',                   ...     % 29    y-component
'sSliceArray\.asSlice\[\d+\]\.sNormal\.dTra',                   ...     % 30    z-component
'sSliceArray.asSlice\[0].dInPlaneRot',                          ...     % 31    The InPlane (InSlice) rotation, so the rotation around the normal vector given by the upper three components
'sGroupArray\.asGroup\[0\]\.dDistFact',                         ...     % 21
'ucUncombImages',                                               ...     % 26 
'sRXSPEC.lGain',                                                ...     % 32 
'sSpecPara\.lPhaseEncodingType',                                ...     % 33    % 1 For Full kSpace Sampling, 2 For Elliptical Weighted Sampling, 3 for Weighted Acquisition
'alTE\[\d+\]',                                                  ...     % 35
'tFree',														...		% 36
'(?<!a)lAverages',												...		% 37	Only search for lAverages excluding alAverages, because this belongs to 'sDiffusion.alAverages.__attribute__.size'
'sWiPMemBlock\.alFree\[(\d){1,2}\]',							...		% 38	All variables set in Special Card + those from above
'sKSpace.ucPhasePartialFourier',								...		% 39
'sKSpace.ucSlicePartialFourier',								...		% 40
'tProtocolName',												...		% 41	How the "Sequence" at the scanner when it was measured was named
'tSequenceFileName',											...		% 42	Which sequence was used
'alTI\[\d+\]',                                                  ...		% 51	TI
'alTR\[\d+\]',                                                  ...		% 52	TR
'sTXSPEC\.asNucleusInfo\[0\]\.flReferenceAmplitude'             ...     % 50
'lTotalScanTimeSec'                                             ...     % 53    Total scan time
'sTXSPEC\.asNucleusInfo\[0\]\.tNucleus'                         ...     % 54    Nucleus
};



% Name the structure entries of ParList like this:
ParList_Assign = { ...
'total_channel_no_measured',                                        ...     % 1     The number of coils with which the data was measured
'total_channel_no_reco',                                            ...     % 27    The number of coils for which the data is reconstructed
'Dwelltimes',                                                       ...     % 18
'LarmorFreq',                                                       ...     % 19
'ThreeD_flag',                                                      ...     % 8
'AsymmetricEcho',                                                   ...     % 11
'InterleavedSliceAcquisition',                                      ...     % 12
'nFreqEnc',                                                         ...     % 2
'nPhasEnc',                                                         ...     % 3
'nPartEnc',                                                         ...     % 6
'nFreqEnc_FinalMatrix',                                             ...     % 4
'nPhasEnc_FinalMatrix',                                             ...     % 5
'nSLC_FinalMatrix',                                                 ...     % 34
'vecSize',                                                          ...     % 9
'RemoveOversampling',                                               ...     % 10
'VoI_Phase',                                                        ...     % 23
'VoI_Read',                                                         ...     % 24
'VoI_Partition',                                                    ...     % 25
'nSLC',                                                             ...     % 7
'SliceThickness',                                                   ...     % 20
'FoV_Phase',                                                        ...     % 13
'FoV_Read',                                                         ...     % 14
'FoV_Partition',                                                    ...     % 22
'PosVOI_Sag',                                                       ...     % 43
'PosVOI_Cor',                                                       ...     % 44
'PosVOI_Tra',                                                       ...     % 45
'SliceNormalVector_VOI_x',                                          ...     % 46
'SliceNormalVector_VOI_y',                                          ...     % 47
'SliceNormalVector_VOI_z',                                          ...     % 48
'InPlaneRotation_VOI',                                              ...     % 49
'Pos_Sag',                                                          ...     % 15
'Pos_Cor',                                                          ...     % 16
'Pos_Tra',                                                          ...     % 17
'SliceNormalVector_x',                                              ...     % 28
'SliceNormalVector_y',                                              ...     % 29
'SliceNormalVector_z',                                              ...     % 30
'InPlaneRotation',                                                  ...     % 31
'SliceGap',                                                         ...     % 21
'SaveUncombined_flag',                                              ...     % 26
'HighGain_flag',                                                    ...     % 32
'Full_ElliptWeighted_Or_Weighted_Acq',                              ...     % 33    1,2 oder 3, ob ihr wirklich...
'TEs',                                                              ...     % 35
'wipMemBlock_tFree',												...		% 36
'nAverages',														...		% 37
'WipMemBlock_alFree',												...		% 38
'PhasePartialFourier',												...		% 39
'SlicePartialFourier',												...		% 40
'tProtocolName',													...		% 41
'tSequenceFileName',												...		% 42
'TIs',                                                              ...		% 51
'TR',                                                               ...		% 52
'RefAmplitude'                                                      ...     % 50
'TotalScanTime'                                                     ...     % 53
'Nucleus'                                                           ...     % 54
};


% Tells function to which format it should convert the found string in the ascconv (remember: all values in the ascconv are strings):
ParList_Convert = { ...
'str2double',                                                       ...     % 1
'char',                                                             ...     % 27
'str2double',                                                       ...     % 18
'str2double',                                                       ...     % 19
'char',                                                             ...     % 8
'char',                                                             ...     % 11
'char',                                                             ...     % 12
'str2double',                                                       ...     % 2
'str2double',                                                       ...     % 3
'str2double',                                                       ...     % 6
'str2double',                                                       ...     % 4
'str2double',                                                       ...     % 5
'str2double',                                                       ...     % 34
'str2double',                                                       ...     % 9
'char',                                                             ...     % 10
'str2double',                                                       ...     % 23
'str2double',                                                       ...     % 24
'str2double',                                                       ...     % 25
'str2double',                                                       ...     % 7
'str2double',                                                       ...     % 20
'str2double',                                                       ...     % 13
'str2double',                                                       ...     % 14
'str2double',                                                       ...     % 22
'str2double',                                                       ...     % 43
'str2double',                                                       ...     % 44
'str2double',                                                       ...     % 45
'str2double',                                                       ...     % 46
'str2double',                                                       ...     % 47
'str2double',                                                       ...     % 48
'str2double',                                                       ...     % 49
'str2double',                                                       ...     % 15
'str2double',                                                       ...     % 16
'str2double',                                                       ...     % 17
'str2double',                                                       ...     % 28
'str2double',                                                       ...     % 29
'str2double',                                                       ...     % 30
'str2double',                                                       ...     % 31
'str2double',                                                       ...     % 21
'char',                                                             ...     % 26
'str2double',                                                       ...     % 32
'str2double',                                                       ...     % 33
'str2double',                                                       ...     % 35
'char',																...		% 36
'str2double',														...		% 37
'str2double',														...		% 38
'char',																...		% 39
'char',																...		% 40
'char',																...		% 41
'char',																...		% 42
'str2double',                                                       ...     % 51
'str2double',                                                       ...     % 52
'str2double'														...		% 50
'str2double'														...		% 53
'char'                                                              ...     % 54
};


% Initialize ParList
for Par_no = 1:numel(ParList_Search)
    eval([ 'ParList.' ParList_Assign{Par_no} ' = ' ParList_Convert{Par_no} '(''0'');' ]);
end

% open file
fid = fopen(file_path,'r');



%% 1. Track down & save element: ASCCONV

begin_found = 0;
ascconv = [];
sLine = 0;

while(sLine > -1)
    
    sLine = fgets(fid); % get string line
    
    if(not(begin_found))                                            % If begin of ascconv not yet found
        
        
        if(not(~contains(sLine,'### ASCCONV BEGIN')))
            begin_found = true;                                     % If current line is begin of ascconv
        else
            continue                                                % If current line is not begin of ascconv --> read next line
        end
        
        
    else                                                            % If begin of ascconv has already been found
        
        if(not(~contains(sLine,'### ASCCONV END ###')))      % If the end was found --> stop while loop
            break
        elseif(not(~contains(sLine,'### ASCCONV BEGIN')))    % In very rare cases there are two 'ASCCONV BEGIN's in the header...
            ascconv = [];
        else
            ascconv = [ascconv; {sLine}];                
            % If current line is not the end --> read in line and save it
        end
        
    end
   
        
end



%% 2. Display error & stop if no Ascconv found

if(not(begin_found))
    display(['Pfui Toifel! You gave me a file without any ascconv, I cannot digest that! Please remember that I am NOT an omnivore.' char(10) ...
             'I will stop here . . .'])
    return
end


%% 3. Convert ascconv

% Convert cell array of strings containing n x 2 entries. The first entries containing the parts before the '=' (pre-=) 
% and the second parts containing the parts after the '=' (post-=)


% Until now ascconv is a cell array of strings (with lets say 348 entries)

% This regexp makes ascconv to a cell array with 348 entries, each of these on its own a cell array of 2 strings
ascconv = regexp(ascconv, '=','split');

% This makes ascconv to a 2*348 = 696x1 cell array of strings; All odd cells contain the parts before the '=', all even cells the part after the '='
ascconv = transpose([ascconv{:}]);

% Now seperate the pre-= and the post-= parts, remove all white spaces before and after the entries.
ascconv = strtrim([ascconv(1:2:end) ascconv(2:2:end)]);

% Now we are happy and have our 348x2 cell array of strings.



%% 4. Search certain entries & Save these

% The following code performs these tasks:
% strfind(...): searches ascconv-ParameterNames (ascconv(:,1)) for the ParList_Search-strings. This results in a cell array, containing [] if in 
% the corresponding cell the Parametername was not found, and [n] if it was found in the corresponding cell on place n of the string;
% not(cellfun(...)): We then search each cell (--> cellfun) if it is empty, and negate the logical output, so that we get the non-empty cells.
% eval(...) We assign the found value to ParList.ParameterName, where ParameterName is determined by ParList_Assign. We also convert the values.

for Par_no = 1:numel(ParList_Search)
    TokenAssignString = regexp(ascconv(:,1),ParList_Search{Par_no},'tokens');    
    Par_Logic = not(cellfun('isempty',TokenAssignString));
	TokenAssignString = TokenAssignString(Par_Logic); 
	while(iscell([TokenAssignString{:}]))
		TokenAssignString = [TokenAssignString{:}];
	end

	if(numel(TokenAssignString) > 0)
        TokenAssignString = cell2mat(strcat(TokenAssignString,'+1,'));
		TokenAssignString = ['([' TokenAssignString(1:end-1) '])'];
	else
        TokenAssignString = '';
	end
	if(find(Par_Logic) > 0)
        eval([ 'ParList.' ParList_Assign{Par_no} TokenAssignString ' = ' ParList_Convert{Par_no} '(ascconv(Par_Logic,2));' ]);
	end
end



%% 5. Change & Correct certain values




% Interpret WipMemBlocks
try
	ParList = InterpretWipMemBlock(ParList); %#ok
catch Error
	fprintf('\nCould not InterpretWipMemBlock, due to following error:\n%s',Error.message);
end




% Convert from 0x1 etc. to logicals

% Remove Oversampling
if(~isempty(regexp(ParList.RemoveOversampling, '1(?!.)','once')))			% (?!.) means not followed by anything (anything = .). E.g. 0x10 will not be matched, but 0x1 will.
    ParList.RemoveOversampling = true;
else 
    ParList.RemoveOversampling = false;
end

% Asymmetric Echo
if(~isempty(regexp(ParList.AsymmetricEcho, '1(?!.)','once')))
    ParList.AsymmetricEcho = true;
else
    ParList.AsymmetricEcho = false;
end

% Interleaved Acquisition: Slice Ordering
if(~isempty(regexp(ParList.InterleavedSliceAcquisition, '1(?!.)','once')))
	ParList.InterleavedSliceAcquisition = false;
elseif(~isempty(regexp(ParList.InterleavedSliceAcquisition, '2(?!.)','once')))
	ParList.InterleavedSliceAcquisition = true;
elseif(~isempty(regexp(ParList.InterleavedSliceAcquisition, '4(?!.)','once')))
	ParList.InterleavedSliceAcquisition = 3;				% Single-Shot	   
end

% 3D_flag
if(~isempty(regexp(ParList.ThreeD_flag, '4(?!.)','once')))
    ParList.ThreeD_flag = true;
else
    ParList.ThreeD_flag = false;
end

% 3D_flag
if(~isempty(regexp(ParList.SaveUncombined_flag, '1(?!.)','once')))
    ParList.SaveUncombined_flag = true;
else
    ParList.SaveUncombined_flag = false;
end

% Phase Partial Fourier
if(~isempty(regexp(ParList.PhasePartialFourier, '1(?!.)','once')))
    ParList.PhasePartialFourier = 4/8;
elseif(~isempty(regexp(ParList.PhasePartialFourier, '2(?!.)','once')))
    ParList.PhasePartialFourier = 5/8;
elseif(~isempty(regexp(ParList.PhasePartialFourier, '4(?!.)','once')))
    ParList.PhasePartialFourier = 6/8;
elseif(~isempty(regexp(ParList.PhasePartialFourier, '8(?!.)','once')))
    ParList.PhasePartialFourier = 7/8;
else
    ParList.PhasePartialFourier = 0;
end

% Slice Partial Fourier
if(~isempty(regexp(ParList.SlicePartialFourier, '1(?!.)','once')))
    ParList.SlicePartialFourier = 4/8;
elseif(~isempty(regexp(ParList.SlicePartialFourier, '2(?!.)','once')))
    ParList.SlicePartialFourier = 5/8;
elseif(~isempty(regexp(ParList.SlicePartialFourier, '4(?!.)','once')))
    ParList.SlicePartialFourier = 6/8;
elseif(~isempty(regexp(ParList.SlicePartialFourier, '8(?!.)','once')))
    ParList.SlicePartialFourier = 7/8;
else
    ParList.SlicePartialFourier = 0;
end

% Corrections


% In the header, there is always the vector size written with oversampling removed. In the .dat-files, the oversampling is never removed. In the IMA files, it is removed, 
% if ParList.RemoveOversampling=true, otherwise not. Thus: .dat-vecsize always has to be multiplied by 2, IMA only in case of RemoveOversampling=false.
% PROBLEM: SPIRAL IS NEVER (?) OVERSAMPLED. FOR NOW: ONLY REMOVE OVERSAMPLING FOR FULLY SAMPLED DATA SET. BUT THIS IS A HACK.
if(numel(strfind(file_path, '.dat')) > 0 ...
	|| (numel(strfind(file_path, '.IMA')) > 0 && ~ParList.RemoveOversampling && ParList.Full_ElliptWeighted_Or_Weighted_Acq ~= 1))
 	ParList.vecSize = 2 * ParList.vecSize;
    if( isfield(ParList,'WipMemBlockInterpretation') && isfield(ParList.WipMemBlockInterpretation,'Rollercoaster') && isfield(ParList.WipMemBlockInterpretation.Rollercoaster,'sNoADCPointsPerCircle') ...
        && ~isnan(ParList.WipMemBlockInterpretation.Rollercoaster.sNoADCPointsPerCircle) )      % if Rollercoaster
        ParList.Dwelltimes = 2 * ParList.Dwelltimes;    % Because we have oversampling enabled, and the scanner interprets that we oversample in the spectral domain, which we dont...
    end
end
if( numel(strfind(file_path, '.IMA') > 0) ) 
	if( isempty(regexpi(ParList.tProtocolName,'spiral')) && isempty(regexpi(ParList.tSequenceFileName,'spiral')))	% PREVIOUSLY HAD: ParList.Full_ElliptWeighted_Or_Weighted_Acq ~= 4 &&. BUT SOMETIMES CSI IS ALSO WEIGHTED!!! % WITH THAT I ASSUME THAT 
	 	ParList.Dwelltimes = 2 * ParList.Dwelltimes;																													% DATASET IS NOT SPIRAL!
	else																																								% THIS IS A HACK!
		fprintf('\n\nWARNING: I DID  N O T  DOUBLE THE DWELLTIMES AS USUAL.')
		fprintf('\nIF YOUR DATASET IS A CONVENTIONAL, FULLY SAMPLED (no elliptical or acquisition weighting) DATASET,\nTHE RESULTS WILL BE WRONG!')
	end
end


% In case of Single-Slice and Multi-slice, the Partitions and the FinalMatrixSizeSlice are always set to 8, which is quite wrong.
if(~ParList.ThreeD_flag)
    ParList.nPartEnc = 1;
    ParList.nSLC_FinalMatrix = ParList.nSLC;
end



% Total channel number correction
ParList.total_channel_no_measured = numel(ParList.total_channel_no_measured);

if(numel(strfind(file_path, '.IMA')) > 0 && ~ParList.SaveUncombined_flag)
    ParList.total_channel_no_reco = 1;               % Dicom files with SaveUncombined_flag = false, are already coil combined
else
    ParList.total_channel_no_reco = ParList.total_channel_no_measured;
end

% Slice Gap
if(ParList.nSLC > 1)
    ParList.SliceGap = ParList.SliceGap .* ParList.SliceThickness(1:ParList.nSLC-1);    % There is one slice gap less than slices.
else
    ParList.SliceGap = 0;
end


% FoV in Partition Direction
ParList.FoV_Partition = sum(ParList.SliceThickness) + sum(ParList.SliceGap);


% Add field gyromagnetic ratio based on nucleus

if(~isempty(regexp(ParList.Nucleus,'1H','ONCE')))
    ParList.GyroMagnRatioOverTwoPi = 42.57747892 * 10^6;
end










% Consistency Checkings

% Check if All Vectors for Multislice Data have same size.
% It could happen that some are initialized by zeros with dimension [1 1],
% and not set because not available in the header,
% but 4 slices were measured, and some other values have dimension [1 4]. 
% This could lead to problems, e.g. the following:
% SliceNormalVector_x = 0; SliceNormalVector_y = [0.12 0.12]; SliceNormalVector_z = [0.9 0.9];

MaxSizeSliceNormalVecs = max(cat(1,numel(ParList.SliceNormalVector_x),numel(ParList.SliceNormalVector_y),numel(ParList.SliceNormalVector_z))); %#ok
%ParList.SliceNormalVector_x = repmat(ParList.SliceNormalVector_x, [1 1+MaxSizeSliceNormalVecs-numel(ParList.SliceNormalVector_x)]);

for Dim = {'x','y','z'}
    eval([ 'ParList.SliceNormalVector_' Dim{1} ' = repmat(ParList.SliceNormalVector_' Dim{1} ', [1+MaxSizeSliceNormalVecs-numel(ParList.SliceNormalVector_' Dim{1} ') 1]);' ]);
end





%% 5. Postparations

fclose(fid);

end











function ParList = InterpretWipMemBlock(ParList)

	if(~isfield(ParList,'wipMemBlock_tFree'))
		ParList.wipMemBlock_tFree = 0;
	end

	if(  isfield(ParList,'WipMemBlock_alFree') && (numel(ParList.WipMemBlock_alFree) > 1 || (numel(ParList.WipMemBlock_alFree)== 1 && ~(ParList.WipMemBlock_alFree == 0)))  )
		
		ParList.WipMemBlockInterpretation.Prescan = CheckWipMemBlockForPrescan(ParList);
		ParList.WipMemBlockInterpretation.OneDCaipi = CheckWipMemBlockForOneDCaipi(ParList);
		ParList.WipMemBlockInterpretation.TwoDCaipi = CheckWipMemBlockForTwoDCaipi(ParList);
		ParList.WipMemBlockInterpretation.Rollercoaster = CheckWipMemBlockForRollercoaster(ParList);
		
	else
		
		ParList.WipMemBlockInterpretation = 0;
		
	end
	


	
	% Prescan
	function PrescanInterpretation = CheckWipMemBlockForPrescan(ParList)
		
		
		% Info about sizes
		% Try to get Prescans Info from wipMemBlock and ONLINE Info from normal ascconv header and mdh.
		% Prescan Info
		if(numel(ParList.WipMemBlock_alFree) > 55)
			try
				%PrescanInterpretation.nReadEnc = ParList.WipMemBlock_alFree(1);
				PrescanInterpretation.PATREFANDIMASCAN.nPhasEnc = ParList.WipMemBlock_alFree(52);
				PrescanInterpretation.PATREFANDIMASCAN.nPartEnc = ParList.WipMemBlock_alFree(53);
				PrescanInterpretation.PATREFANDIMASCAN.nSLC = ParList.WipMemBlock_alFree(54);
				PrescanInterpretation.PATREFANDIMASCAN.nAverages = ParList.WipMemBlock_alFree(55);
				PrescanInterpretation.NOISEADJSCAN.nReadEnc = ParList.WipMemBlock_alFree(56);
			catch
				fprintf(['\nIs there something different than infos about the Prescans in wipMemBlock.alFree[50-55] ?\n' ...
				'Consider writing the Prescans Info into that array for faster read-in.'])
				ParList.WipMemBlock_alFree = 0;
			end
		end


		% Dwelltimes
		if(numel(ParList.WipMemBlock_alFree) > 58)
			PrescanInterpretation.NOISEADJSCAN.Dwelltime = ParList.WipMemBlock_alFree(57);
			PrescanInterpretation.PATREFANDIMASCAN.Dwelltime = ParList.WipMemBlock_alFree(58);
			PrescanInterpretation.ONLINE.Dwelltime = ParList.WipMemBlock_alFree(59);
		end
		
		
		if(~exist('PrescanInterpretation','var'))
			PrescanInterpretation = 0;
		end
		
		
	end


	function RollercoasterInterpretation = CheckWipMemBlockForRollercoaster(ParList)
		if(numel(ParList.WipMemBlock_alFree) > 59)
			try
				RollercoasterInterpretation.sNoADCPointsPerCircle = ParList.WipMemBlock_alFree(60);
                RollercoasterInterpretation.NoOfTempInt = ParList.WipMemBlock_alFree(1);
                RollercoasterInterpretation.ADCDwellTime = ParList.WipMemBlock_alFree(59);
                RollercoasterInterpretation.SpecDwellTime = RollercoasterInterpretation.ADCDwellTime * RollercoasterInterpretation.sNoADCPointsPerCircle / ...
                                                            RollercoasterInterpretation.NoOfTempInt; % This is just an additional check. We have the info anyway in ParList.Dwelltimes
                
                
			catch errie
				fprintf('\nError: Could not interpret ParList.WipMemBlock_alFree[59] - ParList.WipMemBlock_alFree[59] as Rollercoaster Info.')
				RollercoasterInterpretation.sNoADCPointsPerCircle = NaN;
            	RollercoasterInterpretation.NoOfTempInt = NaN;
                RollercoasterInterpretation.ADCDwellTime = NaN;
                RollercoasterInterpretation.SpecDwellTime = NaN;

			end
		else
			RollercoasterInterpretation.sNoADCPointsPerCircle = NaN;
			RollercoasterInterpretation.NoOfTempInt = NaN;
			RollercoasterInterpretation.ADCDwellTime = NaN;
			RollercoasterInterpretation.SpecDwellTime = NaN;			
		end
		
	end
	
	
	% 1D Caipi
	function OneDCaipiInterpretation = CheckWipMemBlockForOneDCaipi(ParList)

		% Check if info is available in wipmemblock
		OneDCaipiInfoAvail = numel(ParList.WipMemBlock_alFree) > 31 && ParList.WipMemBlock_alFree(31) < 9999 ...
		&& ~(numel(ParList.WipMemBlock_alFree) > 39 && ParList.WipMemBlock_alFree(40) == -1) && sum(ParList.WipMemBlock_alFree(31:40));
	
		if(~OneDCaipiInfoAvail)
			OneDCaipiInterpretation = 0;
			return;
		end

		
		% Find out version
		VerNumber = regexpi(regexpi(ParList.tSequenceFileName,'bs_gh_ghsi_2plus1DCaip_v\d+_\d+_?\d*','match'),'\d+_\d+_?\d*','match');
		if(isempty(VerNumber))
			NewVer = 0;
		else
			VerNumber = VerNumber{:};
			Underscores = strfind(VerNumber,'_'); Fromm = Underscores{1}(1)+1; if(numel(Underscores{1} < 2)); To = Fromm+1; else To = Underscores{1}(2)-1; end
			FeatureNumber = VerNumber{1}(Fromm:To); MainPatchNumber = VerNumber{1}(1:Fromm-2);
			NewVer = str2double(FeatureNumber) > 32 || str2double(MainPatchNumber) > 1; clear Underscores VerNumber FeatureNumber
		end
		
		
		StopLoop = false;
		wipNo = 31;
		SliceAliasingIDs = [];
		FoVShifts_x = [];
		FoVShifts_y = [];

		% Loop always over three consecutive wipMemblock Entries.
		while(~StopLoop && ~(wipNo+2 > 41) && ~(((wipNo - 31)/3 + 1)*4 > ParList.nSLC)) 
			StopLoop = (ParList.WipMemBlock_alFree(wipNo) == -1 || ParList.WipMemBlock_alFree(wipNo+1) == -1 || ParList.WipMemBlock_alFree(wipNo+2) == -1);
			if(StopLoop)
				break;
			end

			SliceAliasingIDs_temp = zeros([1 4]);
			FoVShifts_x_temp = zeros([1 4]);
			FoVShifts_y_temp = zeros([1 4]);

			NumTemp1 = ParList.WipMemBlock_alFree(wipNo);
			NumTemp2 = ParList.WipMemBlock_alFree(wipNo+1);
			NumTemp3 = ParList.WipMemBlock_alFree(wipNo+2);
			

			% The info for the Slicealiasing is given as 1000*SliceID1 + 100*SliceID2 + 10 * SliceID3 + 1 * SliceID4
			% The info about the FoV shifts is given by 10000*FoVShift[SliceID1] + 1000*FoVShift[SliceID2] + 100*FoVShift[SliceID3] + 10*FoVShift[SliceID4]		(OLD VERSION)
			% The info about the FoV shifts is given by 100000000*FoVShift[SliceID1] + 1000000*FoVShift[SliceID2] + 10000*FoVShift[SliceID3] + 100*FoVShift[SliceID4]		(NEW VERSION)

			if(NewVer)
				for BlockNo = 3:-1:0
					SliceAliasingIDs_temp(4-BlockNo) = floor(NumTemp1/10^BlockNo);
					FoVShifts_x_temp(4-BlockNo) = floor(NumTemp2/10^(2*BlockNo));
					FoVShifts_y_temp(4-BlockNo) = floor(NumTemp3/10^(2*BlockNo));

					NumTemp1 = NumTemp1 - floor(NumTemp1/10^BlockNo) * 10^BlockNo;
					NumTemp2 = NumTemp2 - floor(NumTemp2/10^(2*BlockNo)) * 10^(2*BlockNo);
					NumTemp3 = NumTemp3 - floor(NumTemp3/10^(2*BlockNo)) * 10^(2*BlockNo);
				end
				FoVShifts_x_temp = FoVShifts_x_temp / 100; FoVShifts_y_temp = FoVShifts_y_temp / 100;
			else
				for BlockNo = 3:-1:0
					SliceAliasingIDs_temp(4-BlockNo) = floor(NumTemp1/10^BlockNo);
					FoVShifts_x_temp(4-BlockNo) = floor(NumTemp2/10^BlockNo);
					FoVShifts_y_temp(4-BlockNo) = floor(NumTemp3/10^BlockNo);

					NumTemp1 = NumTemp1 - floor(NumTemp1/10^BlockNo) * 10^BlockNo;
					NumTemp2 = NumTemp2 - floor(NumTemp2/10^BlockNo) * 10^BlockNo;
					NumTemp3 = NumTemp3 - floor(NumTemp3/10^BlockNo) * 10^BlockNo;
				end
				FoVShifts_x_temp = FoVShifts_x_temp / 10; FoVShifts_y_temp = FoVShifts_y_temp / 10;
			end

			SliceAliasingIDs = cat(2,SliceAliasingIDs,SliceAliasingIDs_temp);
			FoVShifts_x = cat(2,FoVShifts_x,FoVShifts_x_temp);
			FoVShifts_y = cat(2,FoVShifts_y,FoVShifts_y_temp);

			wipNo = wipNo + 3;

			OneDCaipiInterpretation.SliceAliasingIDs = SliceAliasingIDs;
			OneDCaipiInterpretation.FoVShifts_x = FoVShifts_x;
			OneDCaipiInterpretation.FoVShifts_y = FoVShifts_y;

		end

		clear SliceAliasingIDs_temp FoVShifts_x_temp FoVShifts_y_temp BlockNo NumTemp1 NumTemp2 NumTemp3;


		if(numel(SliceAliasingIDs) > 0 && (numel(unique(SliceAliasingIDs)) < numel(SliceAliasingIDs)))
			OneDCaipiInterpretation.SliceParallelImaging_flag = true;
		else
			OneDCaipiInterpretation.SliceParallelImaging_flag = false;
		end


		if(isfield(OneDCaipiInterpretation,'SliceAliasingIDs'))
			% Convert SliceAliasingIDs to SliceAliasingPattern ([1 2 1 2] --> [1 3; 2 4]);
			SliceAliasingPattern = [];
			ProcessedIDs = zeros(size(OneDCaipiInterpretation.SliceAliasingIDs));

			for CurRep = 1:numel(OneDCaipiInterpretation.SliceAliasingIDs)
				CurID = SliceAliasingIDs(CurRep);

				if(sum(ProcessedIDs == CurID) > 0)
					continue;
				end
				ProcessedIDs(CurRep) = CurID;

				% This only works if all SliceAliasingGroups have the same number of aliased slices (e.g. SliceAliasingIDs = [1 2 1 1] would not be allowed)
				SliceAliasingPattern = cat(1,SliceAliasingPattern,find(OneDCaipiInterpretation.SliceAliasingIDs == CurID));
			end

			OneDCaipiInterpretation.SliceAliasingPattern = SliceAliasingPattern;
			OneDCaipiInterpretation.NoOfMeasSlices = size(OneDCaipiInterpretation.SliceAliasingPattern,1);
			
		end


	end


	




	% 2D Caipi
	function TwoDCaipiInterpretation = CheckWipMemBlockForTwoDCaipi(ParList)
		
		
		TwoDCaipInfoAvail = numel(ParList.WipMemBlock_alFree) > 41 && ParList.WipMemBlock_alFree(41) < 10 && ParList.WipMemBlock_alFree(42) < 10 ...
			&& sum(ParList.WipMemBlock_alFree(41:50));
		TwoDCaipInfoDefective = TwoDCaipInfoAvail && numel(ParList.WipMemBlock_alFree) > 49 && ParList.WipMemBlock_alFree(50) == -1 ...
            && ParList.WipMemBlock_alFree(49) ~= 0;
		
		if(TwoDCaipInfoDefective)
			TwoDCaipiInterpretation = -1;
			fprintf('\nWARNING: 2D-CAIPI seems to have been performed, but I could not read the Pattern.')
			return;
		end
		

		if(TwoDCaipInfoAvail && ParList.WipMemBlock_alFree(41) > 0 && ParList.WipMemBlock_alFree(42) > 0)
			TwoDCaipiInterpretation.Skip_Matrix = zeros([ParList.WipMemBlock_alFree(42) ParList.WipMemBlock_alFree(41)]);
			TwoDCaipiInterpretation.VD_Radius = ParList.WipMemBlock_alFree(43);
			MaxAccess = 50; MaxAccess(MaxAccess > numel(ParList.WipMemBlock_alFree)) = numel(ParList.WipMemBlock_alFree);
			Access = ParList.WipMemBlock_alFree(44:MaxAccess); 
			Access(find(Access == 0,1,'first'):end) = [];
			TwoDCaipiInterpretation.Skip_Matrix(Access) = 1;

			if(numel(TwoDCaipiInterpretation.Skip_Matrix) > 1)
				TwoDCaipiInterpretation.TwoDCaipParallelImaging_flag = true;
			else
				TwoDCaipiInterpretation.TwoDCaipParallelImaging_flag = false;
			end

			clear MaxAccess Access;
		end
		
		if(~exist('TwoDCaipiInterpretation','var'))
			TwoDCaipiInterpretation = 0;
		end
		
		
	end



end

