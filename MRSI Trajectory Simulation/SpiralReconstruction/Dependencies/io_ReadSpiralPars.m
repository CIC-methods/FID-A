function [Out] = io_ReadSpiralPars(SpiralFile,TrajFile,Settings,Out)
%
% read_csi_dat Read in raw data from Siemens
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function can read in MRS(I) data in the Siemens raw file format ".DAT" and performs
% some easy Postprocessing steps like zerofilling, Hadamard decoding, Noise Decorrelation etc.
%
%
% [kSpace, Info] = read_csi_dat(file, DesiredSize,ReadInDataSets)
%
% Input: 
% -         file                    ...     Path of MRS(I) file.
% -         DesiredSize             ...     If you want to perform zerofilling or use only a part of the kspace, set
%                                           DesiredSize to your wanted [kx,ky,kz], e.g. [64 64 1].
% -         ReadInDataSets          ...     
%
% Output:
% -         kSpace                      ...     Output data in k-space. In case of SVS this is zero. size: channel x ROW x COL x SLC x Samples x Averages
% -         Info                        ...     
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: This function uses FFTs to get from k- to image-space. This is mathematically wrong, but Siemens seems to do the same when
% creating DICOMS. The only difference is that the images are flipped left/right and up/down.

% This function expects the input to be of form
% [nCha, nAngInt 


%% 0. Preparations


if(~exist('Out','var'))
    Out = struct;
end
if(~exist('Settings','var'))
   Settings = struct; 
end
if(~isfield(Settings,'IncludeRewinder_flag'))
    Settings.IncludeRewinder_flag = false;
end

Out = supp_UpdateRecoSteps(Out,Settings);


%% Read Parameters


% Read Parameters of D2
Out.Par = read_ascconv(SpiralFile);
% Find out IceProgramPara to properly reshape data
Out.Par.IceProgramPara = read_IceProgramParam(SpiralFile);

% Spiral Specific Parameters
run(TrajFile)
NumOfGradPtsForAllTIs = Out.Par.Dwelltimes/10^3 / 10;    % Dwelltime in us / GRAD_RASTER_TIME in us = Number Of Pts per TI
Out.Par.nTempInt = round(numel(dGradientValues{1}) / NumOfGradPtsForAllTIs);
Out.Par.nAngInts = Out.Par.IceProgramPara.Values(32)/Out.Par.nTempInt;
Out.Par.ADC_dt = Out.Par.IceProgramPara.Values(8);
Out.Par.ADC_OverSamp = 10000/Out.Par.ADC_dt;
Out.Par.TrajPts = NumberOfLoopPoints*Out.Par.ADC_OverSamp;
Out.Par.RewPts = NumberOfBrakeRunPoints*Out.Par.ADC_OverSamp;
Out.Par.TrajTotPts = Out.Par.TrajPts + Out.Par.RewPts;
if(Settings.IncludeRewinder_flag)
    Out.Par.TrajPts = Out.Par.TrajTotPts;
    Out.Par.RewPts = 0;
end

% Find exact vecSize number similar as in sequence
nperiods = 320000 / (10/Out.Par.ADC_OverSamp*Out.Par.TrajTotPts);
if(nperiods > 60)
    concat_periods = floor(nperiods / 60 +1);
    nperiods = floor(nperiods / concat_periods);
end
Out.Par.vecSize = nperiods * concat_periods*Out.Par.nTempInt;




