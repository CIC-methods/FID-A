%% Housekeeping

clearvars
close all;
addpath('./Dependencies')

%% Controls, Files & Inputs

Ctrl.ShowFID_flag = true;
Ctrl.ShowSpec_flag = true;
Ctrl.ShowVecPts_flag = true;
Ctrl.ShowMetMaps_flag = true;

InPath = './Data';

SpiralReco.In.Reco1.Data_file  = [InPath '/meas_MID331_bs_csi_spiral_ase_SPICE_v0_4Test08_D1_20x20_FA90_Avg24_FID66623.dat'];
SpiralReco.In.Reco1.Traj_file  = [InPath '/TrajectoryFiles/20x20_SBW1190_nAI7_nTI1_Slew130p5/SpiralTrajfile_SpatInterleaf0.m'];

    

%% Reco Settings & Perform Reco

Settings.Debug.ShowTrajs = true;
Settings.io_ReadSpiralPars.IncludeRewinder_flag = false;% We could use the rewinders also for reco. However, probably would need to adapt the density compensation?!
Settings.ReadInTraj.GradDelay_x_us = 8;                 % Assume that the gradients are delayed by that amount 
Settings.ReadInTraj.GradDelay_y_us = 8;
Settings.CalcOutTraj.fov_overgrid = 1;                  % Comes from the gridding operator. With the current reco (discrete non-uniform FT) doesnt make sense. 
Settings.CalcOutTraj.OverwriteDataSize_woOvergrid = []; % The Fourier Transform allows to reconstruct to any image-domain grid. You can choose any size.
Settings.NonCartReco.Phaseroll_flag = true;             % k-Space data are not acquired simultaneously, -> each spectral peak acquires phase along traj. Correct that.* 
Settings.NonCartReco.DensComp_flag = true;              % Perform density compensation, because constand-density spirals dont have constant density 
Settings.NonCartReco.DensComp.AutoScale_flag = true;    % This option automatically scales the density compensation function. Only rescales output data. 
Settings.NonCartReco.DensComp.Method = 'SpiralHoge1997AbruptChanges'; % This option specifies which method to calculate density compensation function should be used. So far, only this one works. 
Settings.NonCartReco.CircularSFTFoV_flag = true;        % For circular trajectories our effective FoV is circular. Dont even try to reconstruct the corners!
Settings.NonCartReco.ConjInBegin_flag = false;          % Data are conjugated before FT. This will flip the image and spectra.
Settings.NonCartReco.ConjAtEnd_flag = true;             % Data are conjugated after FT. This will only flip the spectra.
Settings.NonCartReco.Correct4SpatialB0_flag = false;    % Correct for spatial B0-effects. Would need a B0-map for that, which we don't have for this example...

[SpiralReco.Out.Reco1] = op_ReadAndRecoBorjanSpiralData(SpiralReco.In.Reco1.Data_file, SpiralReco.In.Reco1.Traj_file, Settings);

%*: Phaseroll: This phase accumulation along the trajectory causes that off-resonance peaks will get blurred (linear phase along spiral = blurring). Different
% off-resonance peaks will get blurred differently. When we try to correct for that, unfortunately our first and last FID points get a bit screwed and blurred,
% but the metabolites get less blurry. In our images here, we only look at the on-resonance water-peak, and therefore we can only see the artifacts of this
% correction, but not the deblurring effect.


%% Plot VecPoints

if(Ctrl.ShowVecPts_flag)
    Offset = 1;
    Rows = 2;
    Slc = 1;
    
    figure; 
    subplot(Rows,2,1), imagesc(abs(SpiralReco.Out.Reco1.Data(:,:,Slc,1))), colorbar, title('MyReco Vec1')
    subplot(Rows,2,2), imagesc(abs(SpiralReco.Out.Reco1.Data(:,:,Slc,2))), colorbar, title('MyReco Vec2')
    subplot(Rows,2,3), imagesc(abs(SpiralReco.Out.Reco1.Data(:,:,Slc,10))), colorbar, title('MyReco Vec10')
end

%% Plot Spectra

if(Ctrl.ShowSpec_flag)
    
    ppm = compute_chemshift_vector_1_2(SpiralReco.Out.Reco1.RecoPar.LarmorFreq,SpiralReco.Out.Reco1.RecoPar.Dwelltimes(1)/10^9,SpiralReco.Out.Reco1.RecoPar.vecSize);
    Voxels = {[12 13 1], [9 12 1], [10 16 1], [14 15 1], [14 11 1], [11 10 1]};

    figure;
    for i = 1:numel(Voxels)
        subplot(3,2,i);
        hold on
        plot(ppm,abs(fftshift(fft(squeeze(SpiralReco.Out.Reco1.Data(Voxels{i}(1),Voxels{i}(2),Voxels{i}(3),:))))),'r'),
        hold off
        title(['Spec x=' num2str(Voxels{i}(1)) ', y=' num2str(Voxels{i}(2)) ', z=' num2str(Voxels{i}(3))])
    end

end

%% FIDs

if(Ctrl.ShowFID_flag)
    
    time = SpiralReco.Out.Reco1.RecoPar.Dwelltimes(1)/10^9 * 0:(SpiralReco.Out.Reco1.RecoPar.vecSize-1);
    Voxels = {[12 13 1], [9 12 1], [10 16 1], [14 15 1], [14 11 1], [11 10 1]};

    figure;
    for i = 1:numel(Voxels)
        subplot(3,2,i);
        hold on
        plot(time,abs(squeeze(SpiralReco.Out.Reco1.Data(Voxels{i}(1),Voxels{i}(2),Voxels{i}(3),:))),'r'),
        hold off
        title(['FID x=' num2str(Voxels{i}(1)) ', y=' num2str(Voxels{i}(2)) ', z=' num2str(Voxels{i}(3))])
    end

end


