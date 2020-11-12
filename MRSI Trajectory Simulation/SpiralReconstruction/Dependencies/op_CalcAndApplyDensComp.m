function [InStruct,AdditionalOut] = op_CalcAndApplyDensComp(InStruct,NUFTOperator,Settings)
%
% op_ReadAndRecoBorjanSpiralData Read and reconstruct data from Borjan Gagoski's Spiral MRSI Sequence
%
% This function was written by Bernhard Strasser, June 2019.
%
%
% The function can read in Spiral MRSI data in the Siemens raw file format ".DAT" and performs
% the reconstruction of the data (Non-Uniform Slow FourierTransform etc.)
%
%
% [kSpace, Info] = read_csi_dat(file, DesiredSize,ReadInDataSets)
%
% Input: 
% -         InStruct            ...  InStruct-struct with fields (only InTraj.GM and RecoPar have to be present)
%                                           RecoSteps:  All the settings of the performed reconstruction steps
%                                           Par:        The original parameters of the read in data
%                                           RecoPar:    The parameters after reconstruction
%                                           Data:       The non-cartesian data data
%                                           NoiseData:  If available, noise with the same scale and size as the Data is produced. Useful for SNR calculations.
%                                           OutTraj:    The Cartesian trajectory which would correspond to the Fourier transform of the Cartesian output image
%                                                       (if we did k-space gridding, those would be the k-space points we would grid to).
%                                           InTraj:     The (spiral) trajectory with which the data were measured.
% -         NUFTOperator    ...  The operator for going from an image to the NonCartesian k-Space data. Used for the AutoScale option.
% -         Settings        ...  Struct with fields
%                                           Method:         String specifying which method for calculating the DCF should be used.
%                                                           So far, only 'SpiralHoge1997AbruptChanges' implemented.
%                                           AutoScale_flag: When true, the NUFTOperator is used to construct and reconstruct an image of ones. The ground truth
%                                                           and reconstruction are compared in scale, and the ratio is applied to the DCF.
%                                                           After that, the DCF should provide properly scaled images.
%
% Output:
% -         InStruct        ...  Same as for input
% -         AdditionalOut   ...  Struct with additional output, e.g. Density Compensation Function, etc.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: ???

% Further remarks:



%% 0. Preparations

if(~exist('Settings','var'))
    Settings = struct;
end

if(~isfield(Settings,'Method'))
   Settings.Method = 'SpiralHoge1997AbruptChanges';    
end
if(~isfield(Settings,'AutoScale_flag'))
   Settings.AutoScale_flag = true;    
end


%% Calculate Density Compensation According to Hoge1997 - Abrupt Changes

if(strcmpi(Settings.Method,'SpiralHoge1997AbruptChanges'))
    v1 = InStruct.InTraj.GM;
    DCFPreG = zeros([size(v1,2) size(v1,3)]);
    for SpirPts = 2:size(v1,2)
        DCFPreG(SpirPts,:) = sqrt( v1(1,SpirPts,:).^2 + v1(2,SpirPts,:).^2 ) .* ...
        abs( sqrt( v1(1,SpirPts,:).^2 + v1(2,SpirPts,:).^2 ) - sqrt( v1(1,SpirPts-1,:).^2 + v1(2,SpirPts-1,:).^2 ) );
    end
    DCFPreG(isnan(DCFPreG)) = 0;
else
   st = dbstack;
   fprintf('\nWarning in %s: Did not recognize method for calculating density compensation function ''%s''.\nUse DCF = 1 instead.\n',st(1).name,Settings.Method)
   DCFPreG = ones(size_MultiDims(InStruct.InTraj.GM,[2 3]));
end
    
%% Scale DCF
    
if(~Settings.AutoScale_flag)
%         FudgeFactor = 1.2743;     % For old trajectory
%         FudgeFactor = 0.00051078;         % For new trajectory
    FudgeFactor = 1.9634;
    Scale = max(DCFPreG(:))*2*InStruct.RecoPar.fov_overgrid^2/FudgeFactor;
% I dont know what these factors are. The 2*SpSpice.SimPar.fov_overgrid^2 I guessed. The FudgeFactor I got by inputting a image of ones
% and seeing how it was scaled...

else
    OnesData = ones(InStruct.RecoPar.DataSize(1:2));
    OutOnesData = abs(NUFTOperator'*(DCFPreG(:) .* (NUFTOperator*OnesData(:)))*size(InStruct.OutTraj.GM(:,:),2));
    OutOnesData(OutOnesData == 0) = NaN;
    Scale = nanmean(OutOnesData);
end
DCFPreG = DCFPreG/Scale;


%% Apply DCF

if(isfield(InStruct,'Data'))
    InStruct.Data = InStruct.Data .* myrepmat(DCFPreG(:),size(InStruct.Data));
end
if(isfield(InStruct,'NoiseData'))
    InStruct.NoiseData = InStruct.NoiseData .* myrepmat(DCFPreG(:),size(InStruct.NoiseData));
end
    

%% Postparations

if(nargout > 1 && exist('DCFPreG','var'))
    AdditionalOut.DCFPreG = DCFPreG;
end

InStruct = supp_UpdateRecoSteps(InStruct,Settings);


