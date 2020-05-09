function [DataStruct] = sim_CalcCartTraj(DataStruct,Settings)
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

if(~exist('Settings','var'))
    Settings = struct;
end
if(~isfield(Settings,'fov_overgrid'))
    Settings.fov_overgrid = 1;
end
if(~isfield(DataStruct,'RecoPar'))
    DataStruct.RecoPar = DataStruct.Par;
end

DataStruct.RecoPar.fov_overgrid = Settings.fov_overgrid;

if(isfield(Settings,'OverwriteDataSize_woOvergrid') && ~isempty(Settings.OverwriteDataSize_woOvergrid))
    DataStruct.RecoPar.DataSize = Settings.OverwriteDataSize_woOvergrid;    
else
    DataStruct.RecoPar.DataSize = [DataStruct.RecoPar.nFreqEnc DataStruct.RecoPar.nPhasEnc DataStruct.RecoPar.nPartEnc*DataStruct.RecoPar.nSLC ...
                                   DataStruct.RecoPar.vecSize DataStruct.RecoPar.total_channel_no_measured];
end
if(DataStruct.RecoPar.fov_overgrid > 1)
    DataStruct.RecoPar.DataSize(1:2) = DataStruct.RecoPar.DataSize(1:2) * DataStruct.RecoPar.fov_overgrid;   % For 3D Trajectories, this would be wrong! 
end                                                                                                          % Would also need to modify 3rd dim!

%% Calculate Cartesian Trajectory

FoV = DataStruct.RecoPar.FoV_Read/1000*DataStruct.RecoPar.fov_overgrid; % in m
DeltaGM = 10^9/(FoV*DataStruct.RecoPar.GyroMagnRatioOverTwoPi);           % in mT/m * us

% Calculate a Grid for PhaseEncoding Steps
[bla_x, bla_y] = find(ones(DataStruct.RecoPar.DataSize(1)));

DataStruct.OutTraj.GM(1,1,:) = bla_x - floor(DataStruct.RecoPar.DataSize(1)/2) - 1; 
DataStruct.OutTraj.GM(2,1,:) = bla_y - floor(DataStruct.RecoPar.DataSize(1)/2) - 1; 

DataStruct.OutTraj.GM = DataStruct.OutTraj.GM * DeltaGM;

DataStruct.OutTraj.maxR = DeltaGM*DataStruct.RecoPar.DataSize(1)/2;


%% Normalize Trajectory

DataStruct.OutTraj.GM = DataStruct.OutTraj.GM/(DataStruct.OutTraj.maxR*2);


%% PostParations

DataStruct = supp_UpdateRecoSteps(DataStruct,Settings);


