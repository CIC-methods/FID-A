function [kSpace] = reshape_BorjanSpiralData(kSpace, mdhInfo, Settings)

%% Prep
tic
fprintf('\n\nReshaping data\t\t\t...')

if(~exist('Settings','var'))
   Settings = struct; 
end

if(~isfield(Settings,'ProduceNoiseThroughAvgs_flag'))
    Settings.ProduceNoiseThroughAvgs_flag = true;
end

kSpace = supp_UpdateRecoSteps(kSpace,Settings);


%% 
if(exist('mdhInfo','var') && ~isempty(mdhInfo))
    kSpace2 = reshape_Siemens_dat(kSpace.Data,mdhInfo,[11 9 2 1 7 4 5],true); %,[1 7 4 5 11 9 2]); 
    kSpace.Data = kSpace2.ONLINE{1}; clear kSpace2
else
    kSpace.Data = permute(kSpace.Data, [11 1 3 2 6 4 5 7 8 9 10]);
end
kSpace.Data = single(kSpace.Data);
clear mdhInfo

% Reshape to [nTempInt x nAngInt x samples x nADCs x nCha x nAvg x nPart x nSlc]
CurSize = size(kSpace.Data); 
CurSize = cat(2,CurSize,ones([1 8-numel(CurSize)]));
kSpace.Data = reshape(kSpace.Data,[kSpace.Par.nTempInt CurSize(1)/kSpace.Par.nTempInt CurSize(2:end)]);
% kSpace.Data = reshape(kSpace.Data,[CurSize(1)/kSpace.Par.nTempInt kSpace.Par.nTempInt CurSize(2:5) prod(CurSize(6:7))]);


% Reshape to [nTempInt x nAngInt x samples*nADCs x nCha x nAvg x nPart x nSlc]
% Cut those ADC-points that are noise only away
ADC_Points = kSpace.Par.vecSize * kSpace.Par.TrajTotPts/kSpace.Par.nTempInt;        % 
Size_D2 = size(kSpace.Data);
kSpace.Data = reshape(kSpace.Data,[Size_D2(1:2) prod(Size_D2(3:4)) Size_D2(5:end)]);


% Cut together the useful data from the different temp interleaves
% (for the first TI the points 1:ADC_Points are useful, for the second TI the points kSpace.Par.TrajTotPts/kSpace.Par.nTempInt+1:ADC_Points+kSpace.Par.TrajTotPts/kSpace.Par.nTempInt are useful,
% for the third TI the points 2*kSpace.Par.TrajTotPts/kSpace.Par.nTempInt+1:ADC_Points+2*kSpace.Par.TrajTotPts/kSpace.Par.nTempInt are useful etc.
PtNumber = size(kSpace.Data,3) - ADC_Points;
LastPts = kSpace.Data(:,:,end-PtNumber+1:end,:,:,:,:);
kSpace.Data = kSpace.Data(:,:,1:ADC_Points,:,:,:,:);
for curTempInt = 1:kSpace.Par.nTempInt
    StartPt = (curTempInt-1)/kSpace.Par.nTempInt*kSpace.Par.TrajTotPts;        
%     EndPt = ADC_Points - StartPt;
    kSpace.Data(curTempInt,:,:,:,:,:,:) = cat(3,kSpace.Data(curTempInt,:,StartPt+1:ADC_Points,:,:,:,:), LastPts(curTempInt,:,1:StartPt,:,:,:,:));
end
clear LastPts

% Reshape to [nTempInt x nAngInt x nTrajPoints x vecSize x nCha x nAvg x nPart x nSlc]
% Cut the rewinder away, and save the final size
kSpace.Data = reshape(kSpace.Data,[Size_D2(1:2) kSpace.Par.TrajTotPts kSpace.Par.vecSize/kSpace.Par.nTempInt Size_D2(5:end)]);
kSpace.Data = kSpace.Data(:,:,1:kSpace.Par.TrajPts,:,:,:,:,:);


% % DEBUG: ONLY USE FIRST FID POINT
% kSpace.Data = kSpace.Data(:,1,:,1,1,1);
% kSpace.Par.vecSize = 1; kSpace.Par.nTempInt = 1; kSpace.Par.total_channel_no_measured = 1; kSpace.Par.vecSize = 1;


% Define Noise
if(~isfield(kSpace,'NoiseData') && Settings.ProduceNoiseThroughAvgs_flag && size(kSpace.Data,6) > 1)
    kSpace.NoiseData = kSpace.Data(:,:,:,:,:,2,:,:) - kSpace.Data(:,:,:,:,:,1,:,:);
    kSpace.NoiseData = kSpace.NoiseData * sqrt(size(kSpace.Data,6)/2); % Scale the noise to be the same as the noise in the actual data is
else
    kSpace.NoiseData = [];
end
    
% Add averages together
kSpace.Data = squeeze_single_dim(sum(kSpace.Data,6),6);  % squeeze causes non-problematic error in case there is no part/slc
% kSpace.Data = kSpace.Data(:,:,:,:,:,1,:,:);

% Permute from [nTempInt x nAngInt x nTrajPoints x vecSize x nCha x nPart x nSlc] 
%           to [nTrajPoints x nAngInt x nPart x nSlc x nTempInt x vecSize x nCha]
kSpace.Data = permute(kSpace.Data,[3 2 6 7 1 4 5]);
if(isfield(kSpace,'NoiseData') && ~isempty(kSpace.NoiseData))
    kSpace.NoiseData = permute(kSpace.NoiseData,[3 2 6 7 1 4 5]);
end

% Reshape to [nTrajPoints x nAngInt x nPart x nSlc x nTempInt*vecSize x nCha]
Size = size(kSpace.Data); Size = cat(2,Size,ones([1 7-numel(Size)]));
kSpace.Data = reshape(kSpace.Data,[Size(1:4) prod(Size(5:6)) Size(7:end)]);
if(isfield(kSpace,'NoiseData') && ~isempty(kSpace.NoiseData))
    kSpace.NoiseData = reshape(kSpace.NoiseData,[Size(1:4) prod(Size(5:6)) Size(7:end)]);
end

kSpace.Par.DataSize = size(kSpace.Data);


%% Postparations

% Output.Data = kSpace.Data;
% if(isfield(kSpace,'NoiseData') && ~isempty(kSpace.NoiseData))
%     Output.NoiseData = kSpace.NoiseData;
% end
% clear kSpace
% Output.kSpace.Par = kSpace.Par;

%% The End

fprintf('\n\t\t\t\t...took\t%10.6f seconds',toc)


