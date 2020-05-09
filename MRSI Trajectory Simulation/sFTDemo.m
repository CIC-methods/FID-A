%% Description
% This script shows the non-Cartesian Fourier transform using a slow Fourier transform algorithm (sFT). This algorithm simply calculates all the terms
% exp(2pi <k_j,r_l>) for all j and l and stores them in a matrix. This matrix can then be used as an operator to go from an image to the non-Cartesian k-space data,
% or by using the adjoint operator to transform non-Cartesian k-space data to the corresponding image. In the latter case a density compensation function must be
% applied before the sFT. The algorithm is then tested for three cases:
% 1) By defining a groundtruth image, transforming that image to non-Cartesian k-space, and then transforming it back to image domain.
% This tests whether the operator (DCF * sFT') is indeed the approximate inverse of the operator sFT.
% 2) By defining a k-space where all values are 0, except of the k-space center. Since we know the result in image-domain, i.e. a uniform image, we can test if
% our algorithm is really doing a Fourier transform and not any other transform fulfilling (DCF * sFT') ~ inv(sFT).
% 3) By defining a Gaußian k-space, where each k-space point has the value exp(-(x^2+y^2)/sigma^2). Since this function is Fourier transformable analytically,
% we know the groundtruth, and can compare our reconstruction to that.
% Written by Bernhard Strasser, 2019-05-23.

%% Housekeeping

clearvars
close all



%% Define Parameters
% Normally these parameters would be read in from the data file
% Define them here, instead

% Input Parameters (Non-Cartesian k-Space Trajectory)
In.Par.GradDelay_x = 0;      % Due to Eddy currents and other effects, the real trajectory will be slightly delayed in comparison to the ideal trajectory
In.Par.GradDelay_y = 0;      % These variables account for that problem. Units: [us]
In.Par.ADC_OverSamp = 2;     % Our gradient time grid is 10 us, but our ADC-sampling points might be 5 us. In that case set ADC_OverSamp = 2
In.Par.nAngInts = 15;        % Number of angular interleaves
In.Par.FoV_Read = 220;       % FoV in mm

% Output Parameters (Cartesian Image)
Out.Par.DataSize = [80, 80, 1, 1, 1];
Out.Par.fov_overgrid = 1;
Out.Par.FoV_Read = 220;
Out.Par.GyroMagnRatioOverTwoPi = 42.57747892 * 10^6;


%% Create GroundTruth Image

Image_GroundTruth = ones(Out.Par.DataSize);
% Image_GroundTruth = zeros(Out.Par.DataSize); Image_GroundTruth(floor(Out.Par.DataSize(1)/2)+1,floor(Out.Par.DataSize(2)/2)+1) = 1;
% Image_GroundTruth = (EllipticalFilter(ones(size(Image_GroundTruth)),[1 2],[1 1 1 4],1));


%% Get Spiral Trajectory

Traj_file = './SpiralTrajfile_SpatInterleaf0.m';
[In.Traj,In.Par] = io_ReadSpiralTraj(In.Par,Traj_file);
[Out.Traj,Out.Par] = sim_CalcCartTraj(Out.Par);
In.Traj.maxR = Out.Traj.maxR;
In.Traj.GM = In.Traj.GM/(In.Traj.maxR*2);
Out.Traj.GM = Out.Traj.GM/(Out.Traj.maxR*2);
par.nAngInts = 25;
%testTraj = Rosette(par);

%% Calculate FourierTransform Operator

SizeData_k = size_MultiDims(In.Traj.GM,[2 3]); 
SizeData_k = cat(2,SizeData_k,ones([1 5-numel(SizeData_k)]));

tic
sft2_Oper = sft2_Operator(transpose(squeeze(Out.Traj.GM(:,:))*Out.Par.DataSize(1)), transpose(In.Traj.GM(:,:)), 1);
%rosetteFTOperator = sft2_Operator(transpose(testTraj(:,:)), transpose(squeeze(Out.Traj.GM(:,:))*Out.Par.DataSize(1)), 0);
toc

% Restrict to circular FoV
FoVMask = EllipticalFilter(ones(Out.Par.DataSize(1:2)),[1 2],[1 1 1 Out.Par.DataSize(1)/2-1],1); 
FoVMask = FoVMask(:);
sft2_Oper(:,~logical(FoVMask(:))) = 0;
rosetteFTOperator(~logical(FoVMask),:) = 0;
clear FoVMask;


%% Calculate Density Compensation Function

% From Hoge1997, abrupt changes
v1 = In.Traj.GM;
DCFPreG = zeros([size(v1,2) size(v1,3)]);
for SpirPts = 2:size(v1,2)
    DCFPreG(SpirPts,:) = sqrt( v1(1,SpirPts,:).^2 + v1(2,SpirPts,:).^2 ) .* ...
    abs( sqrt( v1(1,SpirPts,:).^2 + v1(2,SpirPts,:).^2 ) - sqrt( v1(1,SpirPts-1,:).^2 + v1(2,SpirPts-1,:).^2 ) );
end
DCFPreG(isnan(DCFPreG)) = 0;
DCFPreG(1,:) = DCFPreG(2,:);    % From the formula from Hoge1997, the first DCF-point would be 0. But that causes a small problem if 
                                % I define my input k-space being all 0's except of the first k-space point
                                
% Rescale DCF
OnesData = ones(Out.Par.DataSize(1:2));
OutOnesData = abs(sft2_Oper'*(DCFPreG(:) .* (sft2_Oper*OnesData(:)))*size(Out.Traj.GM(:,:),2));
OutOnesData(OutOnesData == 0) = NaN;
Scale = nanmean(OutOnesData);

DCFPreG = DCFPreG/Scale;

clear OnesData OutOnesData Scale SpirPts



%% TEST1: Construct spiral k-Space data, and then reconstruct it
% The operator we calculated does exactly that: Going from iSpace --> kSpace
In.Data = sft2_Oper*Image_GroundTruth(:);
% PSF = reshape(In.Data, Out.Par.DataSize);

kSpacePoints = size(testTraj,2)*size(testTraj,3);
%PSF = rosetteFTOperator*ones(kSpacePoints, 1);
%PSF = reshape(PSF, Out.Par.DataSize);
[xq,yq] = meshgrid(1:.5:80, 1:.5:80);

figure;
vq = griddata(1:80, 1:80, real(PSF), xq, yq);
subplot(1,2,1), mesh(xq, yq, vq), axis square, title('Point Spread Function')
subplot (1,2,2), plot(testTraj(1,:,1), testTraj(2,:,1)), axis square, title("k space trajectory"), hold on
for i = 2:par.nAngInts
    plot(testTraj(1,:,i), testTraj(2,:,i))
end
hold off;

figure;
PSF = sft2_Oper'*ones(6810,1);
PSF = reshape(PSF, Out.Par.DataSize);
vq = griddata(1:80, 1:80, real(PSF), xq, yq);
subplot(1,2,1), mesh(xq, yq, vq), axis square, title("Point Spread Function")
subplot(1,2,2), plot(In.Traj.GM(1,:,1), In.Traj.GM(2,:,1)), axis square, title("k space trajectory"), hold on
for i = 2:In.Par.nAngInts
    plot(In.Traj.GM(1,:,i), In.Traj.GM(2,:,i))
end

% Reconstruct image data from spiral k-space data
% We cannot invert our operator (the one that does iSpace --> kSpace), this is an ill-posed problem. 
% Instead, approximate the inverse operator by density compensation and the adjoint operator.

Out.Data = (sft2_Oper'*(In.Data .* DCFPreG(:))) * size(Out.Traj.GM(:,:),2);

Out.Data = reshape(Out.Data,Out.Par.DataSize);



%% TEST2: Reconstruct a k-space where only k-space center is nonzero
% Now we know that sft2_Oper'(sft2_Oper(x).*DCF) ~ x, so sft2_Oper'*DCF is an approximate inverse of sft2_Oper. But this still doesnt tell us that
% it is really a non-Cartesian Fourier transform. Let's now transform a spiral k-space to i-space, where only the k-space center is nonzero. We know that this
% should give a uniform image.
% It seems there is a scaling issue here. The result is uniform, but scaled wrongly...


In2 = In;
In2.Data = zeros(size(DCFPreG));
In2.Data(1,:) = 1;
In2.Data = reshape(In2.Data,size(In.Data));

% Reconstruct Test2-Data
Out2 = Out;
Out2.GroundTruthData = ones(Out.Par.DataSize);
Out2.Data = (sft2_Oper'*(In2.Data .* DCFPreG(:))) * size(Out2.Traj.GM(:,:),2);

Out2.Data = reshape(Out2.Data,Out.Par.DataSize);




%% TEST3: Fourier transform Gauß function
% Reconstruct k-space data where we know analytically the Fouriertransform, and compare analytical vs reco

% Define k-space data
DistTokCenter = squeeze(sqrt(In.Traj.GM(1,:,:) .^2 + In.Traj.GM(2,:,:) .^2));
In3 = In;
Sigma = 1/30;
In3.Data = exp(-(DistTokCenter.^2)/(2*Sigma^2))/(2*pi*Sigma^2);
In3.Data = In3.Data(:);

% Define GroundTruth
Out3 = Out;
DistToCentreVoxel = sqrt(squeeze(Out3.Traj.GM(1,1,:).^2 + Out3.Traj.GM(2,1,:).^2)*size(Out3.Traj.GM(:,:),2));
Out3.GroundTruthData = exp(-2*pi^2*Sigma^2*DistToCentreVoxel.^2);
Out3.GroundTruthData = reshape(Out3.GroundTruthData,Out3.Par.DataSize);

% Reconstruct Data
Out3.Data = (sft2_Oper'*(In3.Data .* DCFPreG(:)));      % WHY DONT I NEED MY NORMAL "* size(Out2.Traj.GM(:,:),2)" HERE?!
Out3.Data = reshape(Out3.Data,Out3.Par.DataSize);


%% Plot GroundTruth vs Reco

% Trajectory
% figure; 
% %scatter(squeeze(Out.Traj.GM(1,1,:)),squeeze(Out.Traj.GM(2,1,:)),'b'), hold on   
% %for AngIntNo = 1:In.Par.nAngInts
%     scatter(squeeze(In.Traj.GM(1,:,AngIntNo)), squeeze(In.Traj.GM(2,:,AngIntNo)),'r')
%     plot(squeeze(In.Traj.GM(1,:,AngIntNo)), squeeze(In.Traj.GM(2,:,AngIntNo)),'r')
% end
% hold off, axis square, title('Spiral Trajectory')

% % Test1
% figure;
% MinMax = [min(real(Out.Data(:))) max(real(Out.Data(:)))];
% subplot(2,2,1), plot(abs(In.Data)), axis square, title('Test1: k-Space Data')
% subplot(2,2,3), imagesc(real(Image_GroundTruth),MinMax), colorbar, axis square, title('Ground Truth Image')
% subplot(2,2,4), imagesc(real(Out.Data),MinMax), colorbar, axis square, title('Reco1')
% 
% 
% % Test2
% figure;
% subplot(2,2,1), plot(abs(In2.Data)), axis square, title('Test2: k-Space Data')
% subplot(2,2,3), imagesc(real(Out2.GroundTruthData)), colorbar, axis square, title('Ground Truth Image')
% subplot(2,2,4), imagesc(real(Out2.Data)), colorbar, axis square, title('Reco2')
% 
% 
% % Test3
% figure;
% subplot(2,2,1), plot(abs(In3.Data)), axis square, title('Test3: k-Space Data')
% subplot(2,2,2), imagesc(abs(Out3.Data-Out3.GroundTruthData)), colorbar, axis square, title('Error')
% subplot(2,2,3), imagesc(real(Out3.GroundTruthData)), colorbar, axis square, title('Ground Truth Image')
% subplot(2,2,4), imagesc(real(Out3.Data)), colorbar, axis square, title('Reco3')

