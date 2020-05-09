function [kSpace, Info] = reshape_Siemens_dat(kSpace, Info,DimensionOrder,quiet_flag)

%% Prep

if(~exist('kSpace','var') || ~exist('Info','var'))
    error('Not enough Input: Need at least kSpace data and the mdhInfo')
end
if(~exist('DimensionOrder','var'))
   DimensionOrder = 1:20;
end
if(~exist('quiet_flag','var'))
    quiet_flag = false;
end
if(~quiet_flag)
    tic
    fprintf('\n\nReshaping data\t\t\t...')
end

TotChannels = Info.General.Ascconv.total_channel_no_measured;
EchoLog = strcmp(Info.mdhEntryNames,'EchoIndex');
EchoInd = find(EchoLog);

% Never use the echo dimension, because always use it as cell, so cannot omit this one!
% DimensionOrder(DimensionOrder == EchoInd) = [];
% DimensionOrder(DimensionOrder > EchoInd) = DimensionOrder(DimensionOrder > EchoInd) - 1;
UseDimensions = unique(setdiff(DimensionOrder,EchoInd));
DontUseDimensions = setdiff(1:size(Info.mdhInfo,1),UseDimensions);
% DontUseDimensions(DontUseDimensions == EchoInd) = [];
DimensionOrder(DimensionOrder == EchoInd) = [];
DimensionOrder = cat(2,DimensionOrder,DontUseDimensions);

%% Initialize Data
MeasSets = unique(Info.EvalInfo);
for CurrentMeasNo = 1:numel(MeasSets)
    CurrentMeasSet = MeasSets{CurrentMeasNo};

    % Find all of those meas sets
    FoundMeasSetIndex = strcmp(Info.EvalInfo,CurrentMeasSet);
    DataSize = max(Info.mdhInfo(:,FoundMeasSetIndex),[],2);
    DataSize(DataSize == 0) = 1;    % So that we have no troubles initializing data
    NoOfEchos = DataSize(EchoLog);
    DataSize(1) = TotChannels;    % Since we only read the mdh of the first channel, this is always
    DataSize(DontUseDimensions) = 1;
%     Info.mdhEntryNames(EchoLog) = [];
    kSpace2.(CurrentMeasSet) = cell([1 NoOfEchos]);
    kSpace2.(CurrentMeasSet)(:) = {zeros(DataSize')};    
end

%% Rearrange Data & Info
CurPt = 0;
for CurADC = 1:size(Info.EvalInfo,2)
    CurrentMeasSet = Info.EvalInfo{CurADC};

    Dum = Info.mdhInfo(:,CurADC);
    CurEcho = Dum(EchoInd);
    Dum(DontUseDimensions) = 1;
    Dum(Dum == 0) = 1;

    kSpace2.(CurrentMeasSet){CurEcho}...
    (:,Dum(2),Dum(3),Dum(4),Dum(5),1,Dum(7),Dum(8),:,Dum(10),Dum(11),Dum(12),Dum(13),Dum(14),Dum(15),Dum(16),Dum(17),Dum(18),Dum(19),Dum(20)) = ...
    kSpace(:,CurPt+1:CurPt+Dum(9));

%   Rearrange Info2.mdhInfo

    CurPt = CurPt + Dum(9);

end

kSpace = kSpace2; clear kSpace2;


%% Permute & Do Sequence Specific Things

% Do Sequence Specific Things


% Permute
for CurEcho = 1:numel(kSpace.(CurrentMeasSet))
    kSpace.(CurrentMeasSet){CurEcho} = permute(kSpace.(CurrentMeasSet){CurEcho},DimensionOrder);
end
Info.mdhUsedEntryNames = Info.mdhEntryNames(DimensionOrder(1:numel(UseDimensions)));
Info.mdhUsedEntryNames = cat(2,{['cell{' Info.mdhEntryNames{6} '}']},Info.mdhUsedEntryNames);


%% The End
if(~quiet_flag)
    fprintf('\n\t\t\t\t...took\t%10.6f seconds',toc)
end


