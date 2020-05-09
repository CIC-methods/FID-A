function [Trajectory,Par] = io_ReadSpiralTraj(Par,TrajFile,IncludeRewinder)
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

if(~isfield(Par,'GradDelay_x'))
    Par.GradDelay_x = 0;
end
if(~isfield(Par,'GradDelay_y'))
    Par.GradDelay_y = 0;
end
if(exist('IncludeRewinder','var') && strcmpi(IncludeRewinder,'IncludeRewinder'))
    IncludeRewinder_flag = true;
else
    IncludeRewinder_flag = false;
end    

%% 

run(TrajFile);
AppendZeros = max([ceil(Par.GradDelay_x/10), ceil(Par.GradDelay_y/10)])+2;   % +1 for better interpolation

% The points in our trajectory are GradRasterTime apart from each other.
% Simulate ADC_dt = GradientRaterTime/2
% New interpolation method: The last point has a problem, and I think it's not ok to append a 0 at beginning. So extrapolate
if(IncludeRewinder_flag)
    NoOfTrajPoints = NumberOfLoopPoints + NumberOfBrakeRunPoints; 
else
    NoOfTrajPoints = NumberOfLoopPoints; 
end
CurTraj = dGradientValues{1}(1:NoOfTrajPoints);    

% Trajectory With Gradient Delays (max 20 us in each direction, steps of 1 us)
% Append 2 zeros (--> 20 us) in beginning of trajectory, interpolate to 1us time-grid,
% take the trajectory points according to the gradient delays, go back to time grid of 10/Par.ADC_OverSamp us

GradDelay_x = Par.GradDelay_x;
GradDelay_y = Par.GradDelay_y;
% Append 2 0's --> max of 20 us gradient delay
CurTraj = [zeros([1, AppendZeros]) CurTraj];
% if(AppendZeros > 0)
%     CurTraj = [dGradientValues{1}(end-AppendZeros+1:end) CurTraj];
% end

%removes the last point and extrapolates 2 more
CurTraj = interp1(1:1:NoOfTrajPoints+AppendZeros-1, CurTraj(1:end-1), 1:1:NoOfTrajPoints+AppendZeros+1, 'spline', 'extrap');

% interpolate to fine time grid of 1 us
CurTraj = interp1(1:1:NoOfTrajPoints+AppendZeros+1, CurTraj, 1:0.1:(NoOfTrajPoints+AppendZeros+2-0.1)); 

% Take the correct points, always remove 10 points at end, because we extrapolated too many time points
realTraj = real(CurTraj((AppendZeros*10+1-GradDelay_x):(end-GradDelay_x-10)));
imagTraj = imag(CurTraj((AppendZeros*10+1-GradDelay_y):(end-GradDelay_y-10)));
CurTraj = complex(realTraj, imagTraj);

% undersample trajectory back to 10/Par.ADC_OverSamp us time-grid
CurTraj = CurTraj(1:10/Par.ADC_OverSamp:end);
blaa = -cumsum(CurTraj*dMaxGradAmpl*10/Par.ADC_OverSamp);              % 5 is the ADC_dwelltime in us, GradientRaterTime = 2*ADC_dwelltime
Trajectory.GM(:,:,1) = [real(blaa); imag(blaa)];
Trajectory.GV(:,:,1) = [real(CurTraj); imag(CurTraj)];



clear dGradientValues NoOfTrajPoints NumberOfBrakeRunPoints

% Rotate spiral trajectory
ExtraRot = 0;
InputTraj = Trajectory;
Trajectory.GM = zeros([2, size(InputTraj.GV,2), Par.nAngInts]); 
Trajectory.GV = Trajectory.GM;
for AngIntNo = 1:Par.nAngInts
    CurRotAngle = (AngIntNo-1)/Par.nAngInts*2*pi + ExtraRot/180*pi;
    RotMat = [cos(CurRotAngle) sin(CurRotAngle); -sin(CurRotAngle) cos(CurRotAngle)];   
    Trajectory.GM(:,:,AngIntNo) = RotMat*InputTraj.GM;
    Trajectory.GV(:,:,AngIntNo) = RotMat*InputTraj.GV;
end



