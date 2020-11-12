function [SpiralOut_xt,AdditionalOut] = op_ReadAndRecoBorjanSpiralData(SpiralDataFile,SpiralTrajectoryFile,Settings)
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
% -         SpiralDataFile          ...  The Siemens twix-data file of the spiral sequence   
% -         SpiralTrajectoryFile    ...  The trajectory file containing information about the k-space spiral trajectory.
% -         Settings                ...  Struct with fields
%                                           Debug: The Debug-settings. Subfields:
%                                                           ShowTrajs:  If you want to plot the spiral and Cartesian trajectories to get a feeling how you measured...
%                                           io_ReadSpiralPars:  Settings for the function "io_ReadSpiralPars", see io_ReadSpiralPars.m
%                                           ReadInTraj:         Settings for reading the trajectory, see io_ReadSpiralTraj.m.
%                                           CalcOutTraj:        Settings for calculating the Cartesian trajectory, see sim_CalcCartTraj.m.
%                                           NonCartReco:        Settings for the non-Cartesian MRSI Reco, see op_ReconstructNonCartMRData.m.
%
% Output:
% -         SpiralOut_xt            ...  Output-struct with fields
%                                           RecoSteps:  All the settings of the performed reconstruction steps
%                                           Par:        The original parameters of the read in data
%                                           RecoPar:    The parameters after reconstruction
%                                           Data:       The reconstructed data
%                                           NoiseData:  If available, noise with the same scale and size as the Data is produced. Useful for SNR calculations.
%                                           OutTraj:    The Cartesian trajectory which would correspond to the Fourier transform of the Cartesian output image
%                                                       (if we did k-space gridding, those would be the k-space points we would grid to).
%                                           InTraj:     The (spiral) trajectory with which the data were measured.
%                                           
% -         AdditionalOut           ...  Struct with additional output, e.g. Fourier Transform Operator, Density Compensation Function, etc.
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

if(~isfield_recursive(Settings,'Debug.ShowTrajs'))
    Settings.Debug.ShowTrajs = true;
end

if(~isfield_recursive(Settings,'io_ReadSpiralPars.IncludeRewinder_flag'))
   Settings.io_ReadSpiralPars.IncludeRewinder_flag = false;    
end

if(~isfield_recursive(Settings,'ReadInTraj.GradDelay_x_us'))
   Settings.ReadInTraj.GradDelay_x_us = 8;    
end
if(~isfield_recursive(Settings,'ReadInTraj.GradDelay_y_us'))
   Settings.ReadInTraj.GradDelay_y_us = 8;    
end

if(~isfield_recursive(Settings,'CalcOutTraj.fov_overgrid'))
   Settings.CalcOutTraj.fov_overgrid = 1;    
end

if(~isfield_recursive(Settings,'CalcOutTraj.OverwriteDataSize_woOvergrid'))
   Settings.CalcOutTraj.OverwriteDataSize_woOvergrid = [];    
end


if(~isfield_recursive(Settings,'NonCartReco.CircularSFTFoV_flag'))
   Settings.NonCartReco.CircularSFTFoV_flag = true;    
end
if(~isfield_recursive(Settings,'NonCartReco.Phaseroll_flag'))
   Settings.NonCartReco.Phaseroll_flag = true;    
end
if(~isfield_recursive(Settings,'NonCartReco.DensComp_flag'))
   Settings.NonCartReco.DensComp_flag = true;    
end
if(~isfield_recursive(Settings,'NonCartReco.DensComp.AutoScale_flag'))
   Settings.NonCartReco.DensComp.DensCompAutoScale_flag = true;    
end
if(~isfield_recursive(Settings,'NonCartReco.DensComp.Method'))
   Settings.NonCartReco.DensComp.Method = 'SpiralHoge1997AbruptChanges';    
end
if(~isfield_recursive(Settings,'NonCartReco.ConjInBegin_flag'))
   Settings.NonCartReco.ConjInBegin_flag = false;    
end
if(~isfield_recursive(Settings,'NonCartReco.ConjAtEnd_flag'))
   Settings.NonCartReco.ConjAtEnd_flag = false;    
end
if(~isfield_recursive(Settings,'NonCartReco.Correct4SpatialB0_flag'))
   Settings.NonCartReco.Correct4SpatialB0_flag = false;    
end



%% Read Spiral Parameters

[SpiralOut_xt] = io_ReadSpiralPars(SpiralDataFile,SpiralTrajectoryFile,Settings.io_ReadSpiralPars);


%% Read Spiral Data

% Initial size: [nAngInt*nTempInt x samples x nADCs x nCha x nAvg x nPart x nSlc]

SpiralOut_xt.Data = mapVBVD(SpiralDataFile);
SpiralOut_xt.Data = SpiralOut_xt.Data.image();
SpiralOut_xt = supp_UpdateRecoSteps(SpiralOut_xt,struct(),'mapVBVD');

SpiralOut_xt.Data = single(SpiralOut_xt.Data);
SpiralOut_xt = reshape_BorjanSpiralData(SpiralOut_xt,[]);



%% Read Spiral Trajectory & Cartesian Trajectory

[SpiralOut_xt] = sim_CalcCartTraj(SpiralOut_xt,Settings.CalcOutTraj);

Settings.ReadInTraj.maxR = SpiralOut_xt.OutTraj.maxR;  % Normalize spiral trajectory to maximum of Cartesian trajectory
[SpiralOut_xt] = io_ReadSpiralTraj(SpiralOut_xt,SpiralTrajectoryFile,Settings.ReadInTraj);




%% Spring cleaning

clearvars -except Settings SpiralOut_xt



%% DEBUG: PLOT Spiral Trajectories

if(Settings.Debug.ShowTrajs)
    figure;
    scatter(squeeze(SpiralOut_xt.OutTraj.GM(1,1,:)),squeeze(SpiralOut_xt.OutTraj.GM(2,1,:)),'b'), hold on   
    for AngIntNo = 1:SpiralOut_xt.Par.nAngInts
        scatter(squeeze(SpiralOut_xt.InTraj.GM(1,:,AngIntNo)), squeeze(SpiralOut_xt.InTraj.GM(2,:,AngIntNo)),'r')
        plot(squeeze(SpiralOut_xt.InTraj.GM(1,:,AngIntNo)), squeeze(SpiralOut_xt.InTraj.GM(2,:,AngIntNo)),'r')
    end
    hold off
end



%% Reconstruct Spiral Data

[SpiralOut_xt,AdditionalOut.RecoOperators] = op_ReconstructNonCartMRData(SpiralOut_xt,[],Settings.NonCartReco);




